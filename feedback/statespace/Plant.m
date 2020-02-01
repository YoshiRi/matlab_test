%% Plant

Xlm(:,i) = Adp*Xlm(:,i-1)+Bdp*U(i-1);
Dist = pfront(i)-Cdp*Xlm(:,i);
dV = vfront(i) - Xlm(2,i);
Dists(:,i) = [Dist;dV];
% Dist



invs = Dist/Z0+Snoise(i)*Dist*0; % 1 / Scale(i) 
disparity = BF/Dist+StereoNoise(i)*0; %m Disp (i)

invs = Dist/Z0+Snoise(i)*Dist; % 1 / Scale(i) 
disparity = BF/Dist+StereoNoise(i); %m Disp (i)

%% next input
%U(i) = [F1 F2]*X(i) - dref;