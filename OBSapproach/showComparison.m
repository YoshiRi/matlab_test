figure(10)
plt=plot([0 t(end)],[lam0 lam0],t,X_prop(1,:),t,X_conv1(1,:),'-.',t,X_conv2(1,:),'--',t,X_conv3(1,:),':')
setfigcolor(plt,'grobp')
title('Magnitude Estimation')
xlabel('time [s]')
ylabel('inv depth [1/m]')
xlim([0 END])
grid on
legend('GT','Prop','Conv1','Conv2','Conv3','Location','best')

figure(11)
plt=plot(t,Z,t,X_prop(2,:),t,X_conv1(2,:),'-.',t,X_conv2(2,:),'--',t,X_conv3(2,:),':')
setfigcolor(plt,'grobp')
title('Depth Estimation')
xlabel('time [s]')
ylabel('depth [m]')
xlim([0 END])
grid on
legend('GT','Prop','Conv1','Conv2','Conv3','Location','best')

figure(12)
plt=plot(t,VZ,t,X_prop(3,:),t,X_conv1(3,:),'-.',t,X_conv2(3,:),'--',t,X_conv3(3,:),':')
setfigcolor(plt,'grobp')
title('Velocity Estimation')
xlabel('time [s]')
ylabel('velocity [m/s]')
xlim([0 END])
ylim([-0.2 0.2])
grid on
legend('GT','Prop','Conv1','Conv2','Conv3','Location','best')

%% Error plot
figure(20)
plt=plot(t,abs(X_prop(1,:)-lam0),t,abs(X_conv1(1,:)-lam0),'-.',t,abs(X_conv2(1,:)-lam0),'--',t,abs(X_conv3(1,:)-lam0),':')
setfigcolor(plt,'robp')
title('Magnitude Estimation Error')
xlabel('time [s]')
ylabel('inv depth [1/m]')
xlim([0 END])
grid on
legend('Prop','Conv1','Conv2','Conv3','Location','best')

figure(21)
plt=plot(t,abs(X_prop(2,:)-Z'),t,abs(X_conv1(2,:)-Z'),'-.',t,abs(X_conv2(2,:)-Z'),'--',t,abs(X_conv3(2,:)-Z'),':')
setfigcolor(plt,'robp')
title('Depth Estimation Error')
xlabel('time [s]')
ylabel('depth [m]')
xlim([0 END])
grid on
legend('Prop','Conv1','Conv2','Conv3','Location','best')

figure(22)
plt=plot(t,abs(X_prop(3,:)-VZ'),t,abs(X_conv1(3,:)-VZ'),'-.',t,abs(X_conv2(3,:)-VZ'),'--',t,abs(X_conv3(3,:)-VZ'),':')
setfigcolor(plt,'robp')
title('Velocity Estimation Error')
xlabel('time [s]')
ylabel('velocity [m/s]')
xlim([0 END])
grid on
legend('Prop','Conv1','Conv2','Conv3','Location','best')

%% Save figures
SaveFigPDF(10,['fig/',rename,'MagEst.PDF'])
SaveFigPDF(11,['fig/',rename,'DepEst.PDF'])
SaveFigPDF(12,['fig/',rename,'VelEst.PDF'])
SaveFigPDF(20,['fig/',rename,'MagErr.PDF'])
SaveFigPDF(21,['fig/',rename,'DepErr.PDF'])
SaveFigPDF(22,['fig/',rename,'VelErr.PDF'])
