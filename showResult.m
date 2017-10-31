%% figure
EZ=BF./mDisp;
EVZ = (EZ(2:len) - EZ(1:len-1) )/ST;
EZ2=Z0./Scale;
EVZ2 = (EZ2(2:len) - EZ2(1:len-1) )/ST;

figure(4)
plot(t,X(1,:).','r',[t(1) t(len)],[1/Z0 1/Z0],'b--')
title('Magnitude Estimation')
xlabel('time [s]')
ylabel('inv depth [1/m]')
grid on
legend('EKF','GroundTruth');

figure(5)
plot(t,X(2,:).','r',t,BF./mDisp,'b--',t,EZ2,'m-.',t,Z,'g-.')
title('Depth Estimation')
xlabel('time [s]')
ylabel('depth [m]')
grid on
legend('EKF','Stereo Only','2D Only','Ground Truth');

figure(6)
plot(t,X(3,:).','r',t(1:len-1),EVZ,'b--',t(1:len-1),EVZ2,'m-.',t,VZ,'g-',t(1:len-1),(EVZ2+EVZ)/2,'y:')
title('Velocity Estimation')
xlabel('time [s]')
ylabel('speed of depth [m/s]')
grid on
legend('EKF','Stereo Only','2D Only','Ground Truth','AVE');

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
