%% figure
EZ=BF./mDisp;
EVZ = (EZ(2:len) - EZ(1:len-1) )/ST;
EZ2=Z0./Scale;
EVZ2 = (EZ2(2:len) - EZ2(1:len-1) )/ST;

hfig=figure(4) 
plt = plot(t,X(1,:).','-',[t(1) t(len)],[1/Z0 1/Z0],'--')
setfigcolor(plt,'rg')
title('Magnitude Estimation')
xlabel('time [s]')
ylabel('inv depth [1/m]')
xlim([0 END/2])
grid on
legend('EKF','GroundTruth');
pfig = pubfig(hfig);
pfig.LegendLoc = 'best';
pfig.FigDim = [15 11];
expfig([rename,'MagEstimation'],'-pdf');


hfig=figure(5)
plt=plot(t,Z,'g-.',t,BF./mDisp,'b--',t,EZ2,'y-.',t,X(2,:).','r')
setfigcolor(plt,'gbyr')
title('Depth Estimation')
xlim([0 END/2])
xlabel('time [s]')
ylabel('depth [m]')
grid on
legend('Ground Truth','Stereo Only','2D Only','EKF');
pfig = pubfig(hfig);
pfig.LegendLoc = 'best';
pfig.FigDim = [15 11];
expfig([rename,'DepthEstimation'],'-pdf');


hfig=figure(6)
% plot(t,VZ,'g-',t(1:len-1),EVZ,'b--',t(1:len-1),EVZ2,'m-.',t(1:len-1),(EVZ2+EVZ)/2,'y:',t,X(3,:).','r')
% xlim([0 END/2])
% title('Velocity Estimation')
% xlabel('time [s]')
% ylabel('velocity [m/s]')
% grid on
% legend('Ground Truth','Stereo Only','2D Only','Average','EKF');
% pfig = pubfig(hfig);
% pfig.LegendLoc = 'best';
% pfig.FigDim = [15 11];
% expfig([rename,'VelocityEstimation'],'-pdf');
plt=plot(t,VZ,'g-',t(1:len-1),EVZ,'b--',t(1:len-1),EVZ2,'y-.',t,X(3,:).','r')
setfigcolor(plt,'gbyr')
xlim([0 END/2])
ylim([-1/2 1/2])
title('Velocity Estimation')
xlabel('time [s]')
ylabel('velocity [m/s]')
grid on
legend('Ground Truth','Stereo Only','2D Only','EKF');
pfig = pubfig(hfig);
pfig.LegendLoc = 'best';
pfig.FigDim = [15 11];
expfig([rename,'VelocityEstimation'],'-pdf');


%%
figure(7)
title('KalmanGain')
plt=plot(t,squeeze(KG(:,1,:)),'-',t,squeeze(KG(:,2,:)),'--')
setfigcolor(plt,'gr')
xlabel('time [s]')
ylabel('velocity [m/s]')
grid on

figure(8)
title('Cov')
plot(t,squeeze(P(1,1,:)),'-',t,squeeze(P(2,2,:)),'--',t,squeeze(P(3,3,:)),'-.')
xlabel('time [s]')
ylabel('Covariance')
grid on

%%
Len = length(t);
if check_cinterval
    hfig=figure(9)
    plt=plot([t(1) t(Len)],[1/Z0 1/Z0],'g-',t,X(1,:).','r',t,Xmax(1,:),'b--',t,Xmin(1,:),'b-.')
    setfigcolor(plt,'grbb')
    title('Magnitude Estimation w/ 95% Confidence Interval ')
    xlabel('time [s]')
    ylabel('inv depth [1/m]')
    grid on
    legend('Ground Truth','EKF','UpperBound','LowerBound');
    pfig = pubfig(hfig);
    pfig.LegendLoc = 'best';
    pfig.FigDim = [15 11];
    expfig([rename,'MagEstimationCI'],'-pdf');

    hfig=figure(10)
    plt=plot(t,Z,'g-',t,X(2,:).','r',t,Xmax(2,:),'b--',t,Xmin(2,:),'b-.')
    setfigcolor(plt,'grbb')
    title('Depth Estimation w/ 95% Confidence Interval')
    xlabel('time [s]')
    ylabel('depth [m]')
    grid on
    legend('Ground Truth','EKF','UpperBound','LowerBound');
    xlim([0 END/2])
    pfig = pubfig(hfig);
    pfig.LegendLoc = 'best';
    pfig.FigDim = [15 11];
    expfig([rename,'DepthEstimationCI'],'-pdf');


    hfig=figure(11)
    plt=plot(t,VZ,'g-',t,X(3,:).','r',t,Xmax(3,:),'b--',t,Xmin(3,:),'b-.')
    setfigcolor(plt,'grbb')
    title('Velocity Estimation w/ 95% Confidence Interval')
    xlabel('time [s]')
    ylabel('velocity [m/s]')
    grid on
    legend('Ground Truth','EKF','UpperBound','LowerBound');
    xlim([0 END/2])
    ylim([-1 1])
    pfig = pubfig(hfig);
    pfig.LegendLoc = 'best';
    pfig.FigDim = [15 11];
    expfig([rename,'VelocityEstimationCI'],'-pdf');

end


%%
    hfig=figure(12)
    plt=plot(t,abs(Z-EZ),'b--',t,abs(Z-EZ2),'y-.',t,abs(Z-X(2,:).'),'r')
    setfigcolor(plt,'byr')
    title('Depth Error')
    xlabel('time [s]')
    ylabel('depth error [m]')
    grid on
    legend('Stereo Only','2D Only','EKF');  
    xlim([0 END/2])
    ylim([0 0.02])
    pfig = pubfig(hfig);
    pfig.LegendLoc = 'best';
    pfig.FigDim = [15 11];
    expfig([rename,'DepthError'],'-pdf');


    hfig=figure(13)
    plt=plot(t(1:len-1),abs(VZ(1:len-1)-EVZ),'b--',t(1:len-1),abs(VZ(1:len-1)-EVZ2),'y-.',t,abs(VZ-X(3,:).'),'r')
    setfigcolor(plt,'byr')
    title('Velocity Error ')
    xlabel('time [s]')
    ylabel('velocity error [m/s]')
    grid on
    legend('Stereo Only','2D Only','EKF');
    xlim([0 END/2])
     ylim([0 0.2])
    pfig = pubfig(hfig);
    pfig.LegendLoc = 'best';
    pfig.FigDim = [15 11];
    expfig([rename,'VelocityError'],'-pdf');


%%
mask = ones(size(Z));
mask(mDisp==INFF) =0;
ZstereoErr = sum(abs(Z-EZ).*mask)/sum(mask);
ZscaleErr = sum(abs(Z-EZ2))/len;
ZekfErr= sum(abs(Z-X(2,:).'))/len;

VZstereoErr = sum(abs(VZ(1:len-1)-EVZ).*mask(1:len-1))/sum(mask(1:len-1));
VZscaleErr = sum(abs(VZ(1:len-1)-EVZ2))/(len-1);
VZekfErr= sum(abs(VZ-X(3,:).'))/len;

hfig=figure(14)
val = [ZstereoErr,ZscaleErr,ZekfErr];
b = bar(val);
set( gca, 'XTickLabel', {'Stereo','Scaling','EKF'} )
title('Average Absolute Error in Depth')
    ylabel('depth error[m]')
    pfig = pubfig(hfig);
    pfig.LegendLoc = 'best';
    pfig.FigDim = [15 11];
    expfig([rename,'ZErrorEval'],'-pdf');
    
    hfig=figure(15)
val = [VZstereoErr,VZscaleErr,VZekfErr];
b = bar(val);
set( gca, 'XTickLabel', {'Stereo','Scaling','EKF'} )
title('Average Absolute Error in Velocity')
    ylabel('velocity error [m/s]')
    pfig = pubfig(hfig);
    pfig.LegendLoc = 'best';
    pfig.FigDim = [15 11];
    expfig([rename,'VZErrorEval'],'-pdf');
