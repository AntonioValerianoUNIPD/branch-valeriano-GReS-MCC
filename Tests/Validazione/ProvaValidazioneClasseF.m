clc
clear
close all

%% Input
Data=load('camclay.dat');
epsx=Data(:,1);
realsigmax=Data(:,2);
realp=Data(:,3);
realq=Data(:,4);

% proprietà del materiale
k=0.002;
lambda=0.003;
M=1.2;
E=100000;
e=0;
nu=0.25;
pc=1.5e5; %kPa

% utilities
uno=[1 1 1 0 0 0]';
W=diag([1 1 1 2 2 2]);

% condizioni iniziali
eps=zeros(6,1);
p1=1.5e6;
sigma=[p1 p1 p1 0 0 0]';

% OCR
p0=1/3*sum(sigma(1:3));
OCR=pc/p0;
    Particle = MCC(e,E,nu,M,lambda,k,pc);

%% Risoluzione

J=0;
for i=2:length(epsx)

    epsilon=[epsx(i)-epsx(i-1) 0 0 0 0 0]';
    [sigma, D] = ReturnMapping(Particle, sigma, epsilon);
    eps(1)=eps(1)+epsx(i)-epsx(i-1);

       % plot variables
        p=1/3*sum(sigma(1:3));
        P(i)=p;
        xi=sigma-p*uno;
        WNorm_xi=sqrt(xi'*W*xi);
        Q(i)=sqrt(3/2)*WNorm_xi;
        PC(i)=Particle.pc;

        sigmaxx(i)=sigma(1);
        sigmayy(i)=sigma(2);
        epsxx(i)=eps(1);
        ITER(i)=i;
        EPSV(i)=sum(eps(1:3));
end

%% Plot Risultati

figure (1) %εx-σx
plot(epsxx(2:end),sigmaxx(2:end),'r',LineWidth=1)
hold on
plot(epsx,realsigmax,'.b')
xlabel('εx [-]')
ylabel('σx [KPa]')
title('Diagramma εx-σx')


figure (2) %εx-σy
plot(epsxx(2:end),sigmayy(2:end),'r')
hold on
xlabel('εx [-]')
ylabel('σy [KPa]')
title('Diagramma εx-σy')


figure (3)
    % yielding evolution
plot(P(2:end),Q(2:end),'r',LineWidth=1)
hold on
xlabel('Volumetric Stress p [KPa]')
ylabel('Deviatoric Stress q [KPa]')
title('p-q Diagram')
plot(realp(2:end),realq(2:end),'b.')
        pcyield=PC(2);
        % yielding surface
        pp=linspace(0,pcyield,250);    
        qq=@(pp) sqrt(-pp.*(pp-pcyield).*M^2);
        plot(pp,qq(pp),'b',LineWidth=1.5)
for j=1:length(ITER)
    if mod(ITER(j),250)==0
        pcyield=PC(j);
        pp=linspace(0,pcyield,251);    
        qq=@(pp) sqrt(-pp.*(pp-pcyield).*M^2);
        plot(pp,qq(pp),'b--',LineWidth=0.5)
    else
    end
end
%     % CSL critical state line
    pp=linspace(0,pcyield/2,250);
    qq=@(pp) pp.*M;
    plot(pp,qq(pp),'--k')
legend('p-q','real p-q','Yielding Surface F=0')

sigmaxerror=abs(sigmaxx(2:end)'-realsigmax(2:end));
perror=abs(P(2:end)'-realp(2:end));
qerror=abs(Q(2:end)'-realq(2:end));

figure(5)
tiledlayout(2,1)
nexttile
plot(ITER(2:end),perror,'r')
hold on
plot(ITER(2:end),sigmaxerror,'b')
plot(ITER(2:end),qerror,'g')
xlabel('Iteration')
ylabel('|RealValue-Result|')
legend('p error','sigmax error','q error')
title('Absolute Error')

relperror=abs(perror./P(2:end)');
relqerror=abs(qerror./Q(2:end)');
relsigmaxerror=abs(sigmaxerror./sigmaxx(2:end)');

nexttile
plot(ITER(2:end),relperror,'r')
hold on
plot(ITER(2:end),relsigmaxerror,'b')
plot(ITER(2:end),relqerror,'g')
xlabel('Iteration')
ylabel('|Relative Error|')
legend('p error','sigmax error','q error')
title('Relative Error')

figure(6)
tiledlayout(2,1)
nexttile
plot(epsx(2:end),perror,'r')
hold on
plot(epsx(2:end),sigmaxerror,'b')
plot(epsx(2:end),qerror,'g')
xlabel('εx [-]')
ylabel('|RealValue-Result|')
legend('p error','sigmax error','q error')
title('Absolute Error')

nexttile
plot(epsx(2:end),relperror,'r')
hold on
plot(epsx(2:end),relsigmaxerror,'b')
plot(epsx(2:end),relqerror,'g')
xlabel('εx [-]')
ylabel('|Relative Error|')
legend('p error','sigmax error','q error')
title('Relative Error')
