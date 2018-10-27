% show pole
%close all
P1 = Pole_prop;
P2 = Pole_conv1;
P3 = Pole_conv2;
P4 = Pole_conv3;

cP1 = log(P1)/ST;
cP2 = log(P2)/ST;
cP3 = log(P3)/ST;
cP4 = log(P4)/ST;

% for i = 1:100

%%
xmin_lim = -20;
ymin_lim = -20;
%%
close all
figure(1001)
i = 5
hold on
plt = plot(real(cP1(:,i)),imag(cP1(:,i)),'x',real(cP2(:,i)),imag(cP2(:,i)),'o',real(cP4(:,i)),imag(cP4(:,i)),'d')
setfigcolor(plt,'robp')
grid on
legend('Prop EKF','Conv EKF','Prop OBS','Location','northeast')
% plot([0 0],[ymin_lim -ymin_lim],'LineStyle','--');
hold off
xlabel('Real');ylabel('Imaginary');
xlim([xmin_lim 5]);ylim([ymin_lim -ymin_lim]);
title(['Poles in S plane in t=',num2str(t(i))])
Square_coloring([xmin_lim 0],[0.9,1,1])

figure(1002)
i = 100
hold on
plt = plot(real(cP1(:,i)),imag(cP1(:,i)),'x',real(cP2(:,i)),imag(cP2(:,i)),'o',real(cP4(:,i)),imag(cP4(:,i)),'d')
setfigcolor(plt,'robp')
grid on
legend('Prop EKF','Conv EKF','Prop OBS','Location','northeast')
% plot([0 0],[ymin_lim -ymin_lim],'LineStyle','--');
hold off
xlabel('Real');ylabel('Imaginary');
xlim([xmin_lim 5]);ylim([ymin_lim -ymin_lim]);
title(['Poles in S plane in t=',num2str(t(i))])
Square_coloring([xmin_lim 0],[0.9,1,1])

figure(1003)
i = 200
hold on
plt = plot(real(cP1(:,i)),imag(cP1(:,i)),'x',real(cP2(:,i)),imag(cP2(:,i)),'o',real(cP4(:,i)),imag(cP4(:,i)),'d')
setfigcolor(plt,'robp')
grid on
legend('Prop EKF','Conv EKF','Prop OBS','Location','northeast')
% plot([0 0],[ymin_lim -ymin_lim],'LineStyle','--');
hold off
xlabel('Real');ylabel('Imaginary');
xlim([xmin_lim 5]);ylim([ymin_lim -ymin_lim]);
title(['Poles in S plane in t=',num2str(t(i))])
Square_coloring([xmin_lim 0],[0.9,1,1])

%%
SaveFigPDF(1001,['fig/',rename,'pole1.PDF'])
SaveFigPDF(1002,['fig/',rename,'pole2.PDF'])
SaveFigPDF(1003,['fig/',rename,'pole3.PDF'])


%%
xini = (X_prop(1,1)-lam0)
figure(555)
plt=plot(t,abs(X_prop(1,:)-lam0),t,abs(X_conv1(1,:)-lam0),'-.',t,abs(X_conv3(1,:)-lam0),'--',t,-xini*exp(-3*2*pi*t),':')
setfigcolor(plt,'robg')
title('Magnitude Estimation Error')
xlabel('time [s]')
ylabel('inv depth [1/m]')
xlim([0 END])
grid on
legend('Prop EKF','Conv EKF','Prop OBS','Nominal Convergence','Location','best')
Patches;
xlim([0 2])

%%
x2ini = (X_prop(2,1)-X(2,1))
figure(666)
plt=plot(t,abs(X_prop(2,:)-X(2,:)),t,abs(X_conv1(2,:)-X(2,:)),'-.',t,abs(X_conv3(2,:)-X(2,:)),'--',t,-x2ini*exp(-3*2*pi*t),':')
setfigcolor(plt,'robg')
title('Magnitude Estimation Error')
xlabel('time [s]')
ylabel('inv depth [1/m]')
xlim([0 END])
grid on
legend('Prop EKF','Conv EKF','Prop OBS','Nominal Convergence','Location','best')
Patches;
xlim([0 2])