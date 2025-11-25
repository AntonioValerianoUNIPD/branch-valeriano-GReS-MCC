classdef ModifiedCamClay < handle
  % ELASTOPLASTIC ISOTROPIC material class

  properties (Access = public)
    % Index Ratio
    e
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
    % flag for tabular input parameters
    isTabular
  end

  methods (Access = public)

     % Class constructor method
     function obj = ModifiedCamClay(fID, matFileName, varargin)
        % Calling the function to set the object properties
        if nargin > 2
           obj.isTabular = true;
           obj.readTabMaterialParameters(fID,matFileName,varargin{1});
        else
           obj.readMaterialParameters(fID, matFileName);
        end
     end

    function [status] = initializeStatus(obj, sigma)
      nptGauss = size(sigma,1);
      status = zeros(nptGauss,2);
    end

    % Material stiffness matrix calculation using the object properties
    function [DAll, sigmaOut, status] = getStiffnessMatrix(obj, sigmaIn, epsilon, dt, status, cellID)
      nptGauss = size(sigmaIn,1);
      % Update stress and constitutive tensore with Nonlinear Elastic Predictor algorithm
      [sigmaOut, DAll] = NLinElasticPredictor(obj, sigmaIn, nptGauss, epsilon,cellID);
    end
    
    function [sigmaOut, D] = NLinElasticPredictor(obj, sigmaIn, nptGauss, epsilon, cellID)
        sigma=-sigmaIn;
        if epsilon==0

            for i = 1:nptGauss

                % Elastic Tensor
                I=eye(6);
                uno=[1 1 1 0 0 0]';
                p=1/3*sum(sigma(i,1:3));
                K=(1+obj.e)/obj.k*p;
                mu=3*K*(1-2*obj.nu)/(2*(1+obj.nu));
                D(:,:,i)=K*(uno*uno')+2*mu*(I-1/3*(uno*uno'));
                sigmaOut(i,1:6)=sigmaIn(i,1:6);
            end
            
        else

            for i = 1:nptGauss

                % Preliminaries
                I=eye(6);
                uno=[1 1 1 0 0 0]';
                W=diag([1 1 1 2 2 2]);
                deltaeps=-epsilon(i,:)';
                r=3*(1-2*obj.nu)/(2*(1+obj.nu));
                ce=uno*uno'+2*r*(I-1/3*(uno*uno'));
                deltaepsv=sum(deltaeps(1:3));
                deltaepsd=deltaeps-1/3*deltaepsv*uno;

                %---------- Elastic Trial ----------%
                deltaepsv_e=deltaepsv;
                deltaeps_e=deltaeps;
                p=1/3*sum(sigma(i,1:3));
                avgK=p/deltaepsv_e*(exp((1+obj.e)/obj.k*deltaepsv_e)-1);
                sigma_trial=sigma(i,:)'+avgK*ce*deltaeps_e;
                ptr=1/3*sum(sigma_trial(1:3)); % p_trial
                xi=sigma_trial-ptr*uno;
                WNorm_xi=sqrt(xi'*W*xi);
                qtr=sqrt(3/2)*WNorm_xi; % q_trial

                %---------- Yielding Function F ----------%
                F=(qtr^2)/(obj.M^2)+ptr*(ptr-obj.pc(cellID,i));

                if F<=0

                    %---------- Elastic Behaviour ----------%
                    sigmaOut(i,1:6)=-sigma_trial;
                    % Elastic Stiffness Matrix
                    K=(1+obj.e)/obj.k*p*exp((1+obj.e)/obj.k*deltaepsv_e);
                    psi=(K-avgK)/(deltaepsv_e);
                    D(:,:,i)=avgK*ce+psi*sqrt(W)*ce*(W*deltaeps_e)*uno'*sqrt(W);

                else

                    %---------- Plastic Behaviour ----------%
                    p=1/3*sum(sigma(i,1:3));
                    xi=sigma(i,:)'-p*uno;
                    WNorm_xi=sqrt(xi'*W*xi);
                    q=sqrt(3/2)*WNorm_xi;
                    p_n=p;
                    xi_n=xi;

                    % Compute Consistency Parameter ΔΦ 
                    [p,q,obj.pc(cellID,i),cpar,res,NRiter,resvec] = SolveCpar(p,q,obj,xi,deltaepsv,deltaepsd,r,cellID,i);
                    
                    % Update Values
                    deltaepsv_p=cpar*(2*p-obj.pc(cellID,i));
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
                    sigmaOut(i,1:6)=p*uno+sqrt(2/3)*q*n;
                    sigmaOut(i,1:6)=-sigmaOut(i,1:6);
                    [D(:,:,i)] = PTangentOperator(p,q,obj,K,avgK,cpar,eta,deltaepsv_e,deltaepsd,mean_mu,n,cellID,i);
                end
            end
        end
    end


    function [p,q,pc,cpar,res,NRiter,resvec] = SolveCpar(p,q,obj,xi_n,deltaepsv,deltaepsd,r,cID,i)

        % InitNewtonRaphson
        W=diag([1 1 1 2 2 2]);
        cpar=1.e-10;
        TOL=1.e-10;
        itmax=15;
        NRiter=0;
        pc=obj.pc(cID,i);
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

            % Update solution values
            delta=delta-d_delta;
            p=delta(1); q=delta(2); pc=delta(3); cpar=delta(4);

            % Convergence
            if NRiter==1
                res=norm(g0)/norm(g0);
            else
                res=norm(g)/norm(g0);
            end
            resvec(NRiter)=res;
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

        J(2,1)=-sqrt(3/2)*eta*dmu_dp*(2*n'*W*deltaepsd-(6*eta*cpar/obj.M^2)*WNorm_d_xi);
        J(2,2)=1;
        J(2,3)=-sqrt(3/2)*eta*dmu_dpc*(2*n'*W*deltaepsd-(6*eta*cpar/obj.M^2)*WNorm_d_xi);
        J(2,4)=-sqrt(3/2)*eta*(2*n'*W*deltaepsd*dmu_dcpar-(6*eta*mean_mu/obj.M^2)*(1+cpar*(2*p-pc)/deltaepsv_e)*WNorm_d_xi);

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

    function [D] = PTangentOperator(p,q,obj,K,avgK,cpar,eta,deltaepsv_e,deltaepsd,mean_mu,n,cID,i)
        
        W=diag([1 1 1 2 2 2]);
        I=eye(6);
        uno=[1 1 1 0 0 0]';

        % Parameters a 
        a=1+2*K*cpar+obj.pc(cID,i)*obj.theta*cpar;
        a1=(1+obj.pc(cID,i)*obj.theta*cpar)/a;
        a2=-(2*p-obj.pc(cID,i))/a;
        a3=2*obj.pc(cID,i)*obj.theta*cpar/a;
        a4=obj.theta*obj.pc(cID,i)/K*(2*p-obj.pc(cID,i))/a;
        a5=sqrt(3/2)*eta;
        a6=-3*q*eta/obj.M^2;

        % Parameters b
        b1=-1+(a1/avgK+2*a1*cpar-a3*cpar)*K;
        b2=(2*p-obj.pc(cID,i))+(a2/avgK+2*a2*cpar-a4*cpar)*K;

        % Parameters c
        c=-2*mean_mu*a6*(2*q/obj.M^2)*(1-b2/deltaepsv_e*(sqrt(2/3)*obj.M^2/(2*q)*n'*W*deltaepsd-cpar))-((2*a2-a4)*p-a2*obj.pc(cID,i))*K;
        c1=c^-1*(-2*mean_mu*a6*b1/deltaepsv_e*(sqrt(2/3)*n'*W*deltaepsd-cpar*2*q/obj.M^2)-((a3-2*a1)*p+a1*obj.pc(cID,i))*K);
        c2=c^-1*(2*mean_mu*a5*(2*q)/obj.M^2);

        % Parameters alfa
        alfa1=2*mean_mu*sqrt(2/3)*a5;
        alfa2=K*(a1+a2*c1)-2/3*mean_mu*sqrt(2/3)*a5;
        alfa3=K*a2*c2;
        alfa4=2*mean_mu*a5*sqrt(2/3)*(b1+b2*c1)/deltaepsv_e;
        alfa5=2*mean_mu*a5*sqrt(2/3)*(b2*c2/deltaepsv_e);
        alfa6=2*mean_mu*a6*sqrt(2/3)*(c1+cpar*(b1+b2*c1)/deltaepsv_e);
        alfa7=-2*mean_mu*a5*sqrt(2/3)*(c2+(cpar*b2*c2)/deltaepsv_e);

        % Plastic Tangent Operator
        C1=alfa1*I;
        C2=alfa2*(uno*uno');
        C3=alfa3*sqrt(W)*(uno*n')*sqrt(W);
        C4=alfa4*sqrt(W)*(deltaepsd*uno')*sqrt(W);
        C5=alfa5*sqrt(W)*(deltaepsd*n')*sqrt(W);
        C6=alfa6*sqrt(W)*(n*uno')*sqrt(W);
        C7=alfa7*sqrt(W)*(n*n')*sqrt(W);
        D=C1+C2+C3+C4+C5+C6+C7;

%         dp_deps=(a1+a2*c1)*K*uno+(a2*c2)*K*n;
%         dxi_deps=2*mean_mu*sqrt(2/3)*a5*(I-1/3*(uno*uno')+(b1+b2*c1)/deltaepsv_e*deltaepsd*uno'+b2*c2/deltaepsv_e*deltaepsd*n')+2*mean_mu*sqrt(2/3)*a6*((c1+cpar*(b1+b2*c1)/deltaepsv_e)*n*uno'-(c2+(cpar*b2*c2)/deltaepsv_e)*n*n');
%         D=uno*dp_deps'+dxi_deps;
    end
  end

 methods (Access = private)
    % Assigning material parameters (check also the Materials class)
    % to object properties
    function readMaterialParameters(obj, fID, matFileName)
      tmpVec = readDataInLine(fID, matFileName, 6);
      %
      obj.e = tmpVec(1);
      obj.nu = tmpVec(2);
      obj.M = tmpVec(3);
      obj.lambda = tmpVec(4);
      obj.k = tmpVec(5);
      % Compute Theta
      obj.theta = (1+obj.e)/(obj.lambda-obj.k);
      obj.pc = tmpVec(6);
    end

    function readTabMaterialParameters(obj,fID,fileName,mesh)
       % read parameters in tabular format 
       PoissonRatioFile = readToken(fID,fileName);
       VoidRatioFile = readToken(fID,fileName);
       CLSslopeFile = readToken(fID,fileName);
       VirginCompIndexFile = readToken(fID,fileName);
       SwellRecompIndexFile = readToken(fID,fileName);
       PreconsolidationStressFile = readToken(fID,fileName);
       
       obj.e = setTabularParams(VoidRatioFile,mesh);
       obj.nu = setTabularParams(PoissonRatioFile,mesh);
       obj.M = setTabularParams(CLSslopeFile,mesh);  
       obj.lambda = setTabularParams(VirginCompIndexFile,mesh);
       obj.k = setTabularParams(SwellRecompIndexFile,mesh);
       % Compute Theta
       obj.theta = (1+obj.e)/(obj.lambda-obj.k);
       obj.pc = setTabularParams(PreconsolidationStressFile,mesh);
    end
  end
end
