
%% system matrix
A = [ 1,0,0;0,1,ST;0,0,1];
B = [0; ST^2/2; ST];
Q = B*B.';
% Measurment Covariance
R1 = STEREO_NOISE_S;
R2 = 0.00005;

%% init
Pinit = diag([100,100,100]);
lam0 = 1/Scale(1)/Z(1);
Xinit = [lam0 Z(1) VZ(1)];

P = zeros(3,3,length(t));
X = zeros(3,length(t));
KG = zeros(3,2,length(t));
P(:,:,1) = Pinit;
X(:,1) = Xinit;


%%
BF_ = BF;
for i=2:length(t)
    % Estimate
    Xhat = A * X(:,i-1);
    Phat = A*P(:,:,i-1)*A.' + Q;
    
    % Switch value
    if mDisp(i)== INFF
        R1_ = INFF*INFF; mD=0;
    else
        R1_ = R1; mD = mDisp(i);        
    end
    %H1 = [0 -BF_/Xhat(2)/Xhat(2) 0];
    %H1 = [0 -mD/Xhat(2) 0];
    Phat3=Phat;
    Xhat3=Xhat;
    H3 = [0,ST*Xhat3(3)/Xhat3(2)/Xhat3(2),-ST/Xhat3(2)];
    Kgain3 = Phat3 * H3.' / (H3*Phat3*H3.'+R2);

    % update 3
    Xnew =  Xhat3 + Kgain3*( dScale(i) - 1 + Xhat3(3)/Xhat3(3)*ST);
    Pnew = (eye(3) - Kgain3*H3)*Phat3;

    
    X(:,i) = Xnew;
    P(:,:,i) = Pnew;
    KG(:,1,i)=Kgain;
    %KG(:,2,i)=Kgain2;
end


%% figure
EZ=BF./mDisp;
EVZ = (EZ(2:len) - EZ(1:len-1) )/ST;
EZ2=Z0./Scale;
EVZ2 = (EZ2(2:len) - EZ2(1:len-1) )/ST;

figure(4)
plot(t,X(1,:).','r',[t(1) t(len)],[1/Z0 1/Z0],'b--')
xlabel('time [s]')
ylabel('inv depth [1/m]')
grid on
legend('EKF','GroundTruth');

figure(5)
plot(t,X(2,:).','r',t,BF./mDisp,'b--',t,EZ2,'m-.',t,Z,'g-.')
xlabel('time [s]')
ylabel('depth [m]')
grid on
legend('EKF','Stereo Only','2D Only','Ground Truth');

figure(6)
plot(t,X(3,:).','r',t(1:len-1),EVZ,'b--',t(1:len-1),EVZ2,'m-.',t,VZ,'g-')
xlabel('time [s]')
ylabel('speed of depth [m/s]')
grid on
legend('EKF','Stereo Only','2D Only','Ground Truth');

%%
figure(7)
title('KalmanGain')
plot(t,squeeze(KG(:,1,:)),'-',t,squeeze(KG(:,2,:)),'--')
xlabel('time [s]')
ylabel('speed of depth [m/s]')
grid on

figure(8)
title('Cov')
plot(t,squeeze(P(1,1,:)),'-',t,squeeze(P(2,2,:)),'--',t,squeeze(P(3,3,:)),'-.')
xlabel('time [s]')
ylabel('Covariance')
grid on
