syms tau
Acp = [0 1 0;0 0 1;0 0 -1/tau];
Bcp = [0;0;1/tau];
Ccp = [1 0 0];
syms Ts
sysdp = ss(Acp,Bcp,Ccp,[],Ts)

%% pole analysis
clear
tau = 0.2;
Acp = [0 1 0;0 0 1;0 0 -1/tau];
Bcp = [0;0;1/tau];
Ccp = [1 0 0];

%
Fa =[-0.025, -0.41, 0];
Fb =[-0.030, -0.30, 0];
Fc =[-0.075, -0.25, 0];


PoleA=eig(Acp+Bcp*Fa)
PoleB=eig(Acp+Bcp*Fb)
PoleC=eig(Acp+Bcp*Fc)

%%
figure(1)
clf
hold on
plot(real(PoleA),imag(PoleA),'o')
plot(real(PoleB),imag(PoleB),'x')
plot(real(PoleC),imag(PoleC),'*')
grid on
xlabel('Real')
ylabel('Imag')
legend('Group A','Group B','Group C','Location','Best')
title('Poles of State Feedback ACC controlers')
SaveFigPDF(1,'PoleABC')