% calc input U
Vref = 0;
ierror = ierror + (X(2,i) - Dref);
% ierror = ierror + (Dist - Dref);

% iverror = ierror + (X(3,i) - Vref);

U(:,i) = -ierror * gfb - kfb*Xlm(:,i);
U(:,i) = (X(2,i) - Dref)*Kp + ierror * Ki;
% U(:,i) = (Dist - Dref)*Kp + ierror * Ki;

% div = 2;
% U(:,i) = ((X(3,i) - Vref)*20000 + (X(2,i) - Dref)*20000 + ierror * 500)/div;

%%
%  Kppii = 500/div;
%  Kppip = 20000/div;
%  Kppip2 = 20000/div;
% % 
% Plants = tf([1],[M B K]);
% PIc = tf([Kppip2 Kppii],[1 0]);
% VPlant = PIc*Plants/(1+tf([1 0],[1])*PIc*Plants);
% sys = minreal(Kppip*VPlant/(1+Kppip*VPlant));
% 
% pole(sys)

