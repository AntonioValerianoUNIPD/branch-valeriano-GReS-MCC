classdef MCC < handle
  % ELASTIC ISOTROPIC material class

  properties (Access = public)
    % Index Ratio
    e
    % Elastic Modulus
    E
    % Poisson's ratio
    nu
    % M (CSL slope)
    M
    % Virgin Compression index
    lambda
    % Swell/Recompression index
    k
    % State variable
    theta
    % Preconsolidation Stress
    pc
  end

  methods (Access = public)

     % Class constructor method
     function obj = MCC(e,E,nu,M,lambda,k,pc)
         obj.e=e;
         obj.E=E;
         obj.nu=nu;
         obj.M=M;
         obj.lambda=lambda;
         obj.k=k;
         obj.theta=(1+obj.e)/(obj.lambda-obj.k);
         obj.pc=pc;
         
  
     end
    
     function [sigmaOut, D] = NlinSolv(obj, sigmaIn, epsilon)
        

                        % Preliminaries
                I=eye(6);
                uno=[1 1 1 0 0 0]';
                W=diag([1 1 1 2 2 2]);
                deltaeps=epsilon;
                r=3*(1-2*obj.nu)/(2*(1+obj.nu));
                ce=uno*uno'+2*r*(I-1/3*(uno*uno'));
                deltaepsv=sum(deltaeps(1:3));
                deltaepsd=deltaeps-1/3*deltaepsv*uno;

                %---------- Elastic Trial ----------%
                deltaepsv_e=deltaepsv;
                deltaeps_e=deltaeps;
                p=1/3*sum(sigmaIn(1:3));
                if deltaeps==0
                    avgK=0;
                else
                    avgK=p/deltaepsv_e*(exp((1+obj.e)/obj.k*deltaepsv_e)-1);
                end
                sigma_trial=sigmaIn+avgK*ce*deltaeps_e;
                ptr=1/3*sum(sigma_trial(1:3)); % p_trial
                xi=sigma_trial-p*uno;
                WNorm_xi=sqrt(xi'*W*xi);
                qtr=sqrt(3/2)*WNorm_xi; % q_trial

                %---------- Yielding Function F ----------%
                F=(qtr^2)/(obj.M^2)+ptr*(ptr-obj.pc);

                if F<=0

                    %---------- Elastic Behaviour ----------%
                    sigmaOut=sigma_trial;
                    % Elastic Stiffness Matrix
                    K=(1+obj.e)/obj.k*p*exp((1+obj.e)/obj.k*deltaepsv_e);
                    if avgK==0
                        psi=0;
                    else
                        psi=(K-avgK)/deltaepsv_e;
                    end
                    D=avgK*ce+psi*ce*deltaeps_e*uno';

                else

                    %---------- Plastic Behaviour ----------%
                    p=1/3*sum(sigmaIn(1:3));
                    xi=sigmaIn-p*uno;
                    WNorm_xi=sqrt(xi'*W*xi);
                    q=sqrt(3/2)*WNorm_xi;
                    p_n=p;
                    xi_n=xi;

                    % Compute Consistency Parameter ΔΦ 
                    [p,q,obj.pc,cpar] = SolveCpar(p,q,obj,xi,deltaepsv,deltaepsd,r);

                    % Update Values
                    deltaepsv_p=cpar*(2*p-obj.pc);
                    deltaepsv_e=deltaepsv-deltaepsv_p;
                    avgK=(p-p_n)/deltaepsv_e;
                    K=(1+obj.e)/obj.k*p;
                    mean_mu=r*avgK;
                    eta=(1+(6*mean_mu*cpar)/obj.M^2)^-1;
                    d_xi=xi_n+2*mean_mu*deltaepsd;
                    xi=eta*d_xi;
                    WNorm_xi=sqrt(xi'*W*xi);
                    n=xi/WNorm_xi;

                    % Plastic Stiffness Matrix
                    [D] = PTangentOperator(p,q,obj,K,avgK,cpar,eta,deltaepsv_e,deltaepsd,mean_mu,n);
                    sigmaOut=p*uno+sqrt(2/3)*q*n;
                end
     end


    function [p,q,pc,cpar] = SolveCpar(p,q,obj,xi_n,deltaepsv,deltaepsd,r)

        % InitNewtonRaphson
        W=diag([1 1 1 2 2 2]);
        cpar=0;
        TOL=1.e-10;
        itmax=15;
        NRiter=0;
        pc=obj.pc;
        pc_n=pc;
        delta=[p q pc cpar]'; 
        p_n=p;
        res=2*TOL; % dummy value

        % NR-Cycle
        while res>TOL && NRiter<itmax

            % Update Values
            NRiter=NRiter+1;
            deltaepsv_e=deltaepsv-cpar*(2*p-pc);
            mean_mu=(p-p_n)/deltaepsv_e*r;
            eta=(1+(6*mean_mu*cpar)/obj.M^2)^-1;
            d_xi=xi_n+2*mean_mu*deltaepsd;
            WNorm_d_xi=sqrt(d_xi'*W*d_xi);
            xi=eta*d_xi;

            % NR-Residual g
            if NRiter==1
                g0(1,1)=p-p_n*exp((1+obj.e)*deltaepsv_e/obj.k);
                g0(2,1)=q-sqrt(3/2)*eta*WNorm_d_xi;
                g0(3,1)=pc-pc_n*exp(obj.theta*cpar*(2*p-pc));
                g0(4,1)=q^2/obj.M^2+p*(p-pc);
            else
                g(1,1)=p-p_n*exp((1+obj.e)*deltaepsv_e/obj.k);
                g(2,1)=q-sqrt(3/2)*eta*WNorm_d_xi;
                g(3,1)=pc-pc_n*exp(obj.theta*cpar*(2*p-pc));
                g(4,1)=q^2/obj.M^2+p*(p-pc);
            end

            % MCC-Jacobian J
            [MCC_J] = MCC_Jacobian(p,p_n,q,pc,pc_n,obj,cpar,eta,xi,WNorm_d_xi,mean_mu,r,deltaepsv_e,deltaepsd);

            % SolveDelta 
            if NRiter==1
                d_delta=MCC_J\g0;
            else
                d_delta=MCC_J\g;
            end

            delta=delta-d_delta;
            p=delta(1); q=delta(2); pc=delta(3); cpar=delta(4);

            % Convergence
            if NRiter==1
                res=norm(g0)/norm(g0);
            else
                res=norm(g)/norm(g0);
            end
        end
    end  

    function [J] = MCC_Jacobian(p,p_n,q,pc,pc_n,obj,cpar,eta,xi,WNorm_d_xi,mean_mu,r,deltaepsv_e,deltaepsd)
    
        W=diag([1 1 1 2 2 2]);
        Kstar=(1+obj.e)/obj.k*p_n*exp((1+obj.e)/obj.k*deltaepsv_e);
        pcstar=pc_n*exp(obj.theta*cpar*(2*p-pc));
        WNorm_xi=sqrt(xi'*W*xi);
        n=xi/WNorm_xi;
        if WNorm_xi <= 1.e-7
            n=zeros(6,1);
        end
        
        % row 1
        J(1,1)=1+2*cpar*Kstar;
        J(1,2)=0;
        J(1,3)=-cpar*Kstar;
        J(1,4)=(2*p-pc)*Kstar;

        % row 2
        dmu_dp=(r+2*cpar*mean_mu)/deltaepsv_e;
        dmu_dpc=-(mean_mu*cpar)/deltaepsv_e;
        dmu_dcpar=(2*p-pc)*mean_mu/deltaepsv_e;

        J(2,1)=-sqrt(3/2)*eta*dmu_dp*(2*n'*deltaepsd-(6*eta*cpar/obj.M^2)*WNorm_d_xi);
        J(2,2)=1;
        J(2,3)=-sqrt(3/2)*eta*dmu_dpc*(2*n'*deltaepsd-(6*eta*cpar/obj.M^2)*WNorm_d_xi);
        J(2,4)=-sqrt(3/2)*eta*(2*n'*deltaepsd*dmu_dcpar-(6*eta*mean_mu/obj.M^2)*(1+cpar*(2*p-pc)/deltaepsv_e)*WNorm_d_xi);

        % row 3
        J(3,1)=-2*obj.theta*cpar*pcstar;
        J(3,2)=0;
        J(3,3)=1+obj.theta*cpar*pcstar;
        J(3,4)=-obj.theta*(2*p-pc)*pcstar;

        % row 4
        J(4,1)=2*p-pc;
        J(4,2)=2*q/obj.M^2;
        J(4,3)=-p;
        J(4,4)=0;
    end

    function [D] = PTangentOperator(p,q,obj,K,avgK,cpar,eta,deltaepsv_e,deltaepsd,mean_mu,n)
        
        W=diag([1 1 1 2 2 2]);
        I=eye(6);
        uno=[1 1 1 0 0 0]';

        % Parameters a 
        a=1+2*K*cpar+obj.pc*obj.theta*cpar;
        a1=(1+obj.pc*obj.theta*cpar)/a;
        a2=-(2*p-obj.pc)/a;
        a3=2*obj.pc*obj.theta*cpar/a;
        a4=obj.theta*obj.pc/K*(2*p-obj.pc)/a;
        a5=sqrt(3/2)*eta;
        a6=-3*q*eta/obj.M^2;

        % Parameters b
        b1=-1+(a1/avgK+2*a1*cpar-a3*cpar)*K;
        b2=(2*p-obj.pc)+(a2/avgK+2*a2*cpar-a4*cpar)*K;

        % Parameters c
        c=-2*mean_mu*a6*(2*q/obj.M^2)*(1-b2/deltaepsv_e*(sqrt(2/3)*obj.M^2/(2*q)*n'*deltaepsd-cpar))-((2*a2-a4)*p-a2*obj.pc)*K;
        c1=c^-1*(-2*mean_mu*a6*b1/deltaepsv_e*(sqrt(2/3)*n'*deltaepsd-cpar*2*q/obj.M^2)-((a3-2*a1)*p+a1*obj.pc)*K);
        c2=c^-1*(2*mean_mu*a5*(2*q)/obj.M^2);

        % Parameters alfa
        alfa1=2*mean_mu*sqrt(2/3)*a5;
        alfa2=K*(a1+a2*c1)-2/3*mean_mu*sqrt(2/3)*a5;
        alfa3=K*a2*c2;
        alfa4=2*mean_mu*a5*sqrt(2/3)*(b1+b2*c1)/deltaepsv_e;
        alfa5=2*mean_mu*a5*sqrt(2/3)*(b2*c2/deltaepsv_e);
        alfa6=2*mean_mu*a6*sqrt(2/3)*(c1+cpar*(b1+b2*c1)/deltaepsv_e);
        alfa7=-2*mean_mu*a5*sqrt(2/3)*(c2+cpar*b2*c2/deltaepsv_e);

        % Plastic Tangent Operator
        C1=alfa1*I;
        C2=alfa2*(uno*uno');
        C3=alfa3*uno*n';
        C4=alfa4*deltaepsd*uno';
        C5=alfa5*deltaepsd*n';
        C6=alfa6*n*uno';
        C7=alfa7*(n*n');

        D=C1+C2+C3+C4+C5+C6+C7;
    end
  end
end


