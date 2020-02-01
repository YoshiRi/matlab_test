%%

fname = '_EKF_cont_wonoise'
fname = '_OBS_cont_wonoise'
fname = '_OBS_cont_wnoise'
fname = '_EKF_cont_wnoise'

% fname = '_EKF_obs_wonoise'

%%
SaveFigPDF(9,['fig\MagnitudeEstimation',fname])
SaveFigPDF(10,['fig\Position',fname])
SaveFigPDF(12,['fig\Velocity',fname])
SaveFigPDF(11,['fig\PositionEstimation',fname])
SaveFigPDF(13,['fig\VelocityEstimation',fname])
SaveFigPDF(14,['fig\DistanceControl',fname])

%%

figure(1000)
plot([t(1) t(end)],[dref dref],':',t,pfront-Xlm_EKF(1,:)','-',t,pfront-Xlm_OBS(1,:)','--')
grid on;
legend('Reference','Tracked with EKF','Tracked with OBS','Location','best')
title('Distance reference and tracked result')
xlabel('time [s]')
ylabel('relative position [m]')

figure(1001)
plot(t,vfront,':',t,Xlm_EKF(2,:),'-',t,Xlm_OBS(2,:),'--')
title('Velocity tracking results')
legend('Reference','Tracked with EKF','Tracked with OBS','Location','best')
grid on;
xlabel('time [s]')
ylabel('velocity [m]')

figure(1002)
plot(t,pfront-Xlm_EKF(1,:)',':',t,X_prop(2,:),'-',t,pfront-Xlm_OBS(1,:)',':',t,X_conv3(2,:),'--')
title('Distance estimation')
legend('Ground Truth in EKF','Estimated with EKF','Ground Truth in OBS','Estimated with OBS','Location','best')
grid on;
xlabel('time [s]')
ylabel('relative position [m]')

figure(1003)
plot(t,vfront-Xlm_EKF(2,:)',':',t,X_prop(3,:),'-',t,vfront-Xlm_OBS(2,:)',':',t,X_conv3(3,:),'--')
legend('Ground Truth in EKF','Estimated with EKF','Ground Truth in OBS','Estimated with OBS','Location','best')
grid on;
title('Velocity estimation')
grid on;
xlabel('time [s]')
ylabel('relative velocity [m/s]')

%%
SaveFigPDF(1000,'fig\DistanceControl_comp_wnoise')
SaveFigPDF(1001,'fig\VelocityControl_comp_wnoise')

