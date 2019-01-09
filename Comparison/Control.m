%% Linear motor
M = 2;
B = 100;
K = (4.3*2*pi)^2*M;
%input:i , output:Z

Alm = [0 1;-B/M -K/M];
Blm = [0 ;1/M];
Clm = [1 0];

Alm_ = [Alm [0;0];-1 0 0];
Blm_ = [Blm; 0];
D_ = [0;0;1];

K = place(Alm_,Blm_,[-1;-1+0.3i;-1-0.3i]*2*pi);


%% control design
Zlm = Xlm(1);
dZlm = Xlm(2);

Xlm = Alm*Xlm + Blm*u;
Dist = Z - Zlm;
dV = VZ - dZlm;