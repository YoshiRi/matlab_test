%% Linear motor
M = 2;
B = 100;
K = (4.3*2*pi)^2*M;
sys = tf([1],[M B K]);

%input:i , output:Z

Alm = [0 1;-B/M -K/M];
Blm = [0 ;1/M];
Clm = [1 0];

Alm_ = [Alm [0;0];-1 0 0];
Blm_ = [Blm; 0];
D_ = [0;0;1];

[Alm_d Blm_d]=c2d(Alm_,Blm_,ST);
control_speed = 0.2;%0.2 Hz
K = place(Alm_d,Blm_d,exp(control_speed*[-1;-1+0.3i;-1-0.3i]*2*pi*ST));
kfb = K(1:2);
gfb = K(3);

%% control jissai
Xlm = zeros(2,length(t));
Dists = zeros(2,length(t));
U = zeros(1,length(t));
% Zlm = Xlm(1);
% dZlm = Xlm(2);
% 
% Xlm = Alm*Xlm + Blm*u;
% Dist = Z - Zlm;
% dV = VZ - dZlm;
% 
Dref = 0.4;
% Dref - Dist

%%
Kp = 2000; Ki = 500; 
csys = tf([Kp Ki],[1 0]);
FB = minreal((sys*csys)/(1+sys*csys));
eig(FB)
max(real(eig(FB)))
slowest_control_pole = max(real(eig(FB)))
