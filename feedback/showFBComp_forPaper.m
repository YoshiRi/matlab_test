
close all;
figure(10)
plt=plot([0 t(end)],[lam0 lam0],':',t,X_prop(1,:),t,X_conv1(1,:),'-.',t,X_conv3(1,:),'--')
setfigcolor(plt,'grobp')
title('Magnitude Estimation')
xlabel('time [s]')
ylabel('inv depth [1/m]')
xlim([0 END])
grid on
legend('Ground Truth','Prop EKF','Conv EKF','Prop OBS','Location','best')
Patches;

figure(11)
plt=plot(t,Dists(:,1),':',t,X_prop(2,:),t,X_conv1(2,:),'-.',t,X_conv3(2,:),'--')
setfigcolor(plt,'grobp')
title('Depth Estimation')
xlabel('time [s]')
ylabel('depth [m]')
xlim([0 END])
grid on
legend('Ground Truth','Prop EKF','Conv EKF','Prop OBS','Location','best')
Patches;

figure(12)
plt=plot(t,Dists(:,2),':',t,X_prop(3,:),t,X_conv1(3,:),'-.',t,X_conv3(3,:),'--')
setfigcolor(plt,'grobp')
title('Velocity Estimation')
xlabel('time [s]')
ylabel('velocity [m/s]')
xlim([0 END])
ylim([-0.2 0.2])
grid on
legend('Ground Truth','Prop EKF','Conv EKF','Prop OBS','Location','best')
Patches;
%% Error plot
figure(20)
plt=plot(t,abs(X_prop(1,:)-lam0),t,abs(X_conv1(1,:)-lam0),'-.',t,abs(X_conv3(1,:)-lam0),'--')
setfigcolor(plt,'robp')
title('Magnitude Estimation Error')
xlabel('time [s]')
ylabel('inv depth [1/m]')
xlim([0 END])
grid on
legend('Prop EKF','Conv EKF','Prop OBS','Location','best')
Patches;

figure(21)
plt=plot(t,abs(X_prop(2,:)-Dists(:,1)'),t,abs(X_conv1(2,:)-Dists(:,1)'),'-.',t,abs(X_conv3(2,:)-Dists(:,1)'),'--')
setfigcolor(plt,'robp')
title('Depth Estimation Error')
xlabel('time [s]')
ylabel('depth [m]')
xlim([0 END])
grid on
legend('Prop EKF','Conv EKF','Prop OBS','Location','best')
ylim([0 0.2])
Patches;

figure(22)
plt=plot(t,abs(X_prop(3,:)-Dists(:,2)'),t,abs(X_conv1(3,:)-Dists(:,2)'),'-.',t,abs(X_conv3(3,:)-Dists(:,2)'),'--')
setfigcolor(plt,'robp')
title('Velocity Estimation Error')
xlabel('time [s]')
ylabel('velocity [m/s]')
xlim([0 END])
grid on
legend('Prop EKF','Conv EKF','Prop OBS','Location','best')
ylim([0 0.2])
Patches;
%% Save figures
if savefig
SaveFigPDF(10,['fig/',rename,'MagEst.PDF'])
SaveFigPDF(11,['fig/',rename,'DepEst.PDF'])
SaveFigPDF(12,['fig/',rename,'VelEst.PDF'])
SaveFigPDF(20,['fig/',rename,'MagErr.PDF'])
SaveFigPDF(21,['fig/',rename,'DepErr.PDF'])
SaveFigPDF(22,['fig/',rename,'VelErr.PDF'])
end