% show pole
close all
P1 = Pole_prop;
P2 = Pole_conv1;
P3 = Pole_conv2;
P4 = Pole_conv3;

% for i = 1:100
i = 3


%%
figure(102)
subplot(2,2,1)
hold on
viscircles([0 0],1,'LineStyle','--');
plot(real(P1(:,i)),imag(P1(:,i)),'o')
hold off
grid on
% legend('1st ','2nd ','3rd','4th','Location','northeast')
legend('Poles','Location','northeast')
title('Poles for Prop1')
xlabel('Real');ylabel('Imaginary');
xlim([-1 1.5]);ylim([-2 2]);
%%
subplot(2,2,2)
hold on
viscircles([0 0],1,'LineStyle','--');
plot(real(P2(:,i)),imag(P2(:,i)),'o')
hold off
grid on
% legend('1st ','2nd ','3rd','4th','Location','northeast')
legend('Poles','Location','northeast')
title('Poles for Conv1')
xlabel('Real');ylabel('Imaginary');
xlim([-1 1.5]);ylim([-2 2]); 
%%
subplot(2,2,3)
hold on
viscircles([0 0],1,'LineStyle','--');
plot(real(P3(:,i)),imag(P3(:,i)),'o')
hold off
grid on
% legend('1st ','2nd ','3rd','4th','Location','northeast')
legend('Poles','Location','northeast')
title('Poles for Conv2')
xlabel('Real');ylabel('Imaginary');
xlim([-1 1.5]);ylim([-2 2]);
%% 
subplot(2,2,4)
hold on
viscircles([0 0],1,'LineStyle','--');
plot(real(P4(:,i)),imag(P4(:,i)),'o')
hold off
grid on
% legend('1st ','2nd ','3rd','4th','Location','northeast')
legend('Poles','Location','northeast')
title('Poles for Conv3')
xlabel('Real');ylabel('Imaginary');
xlim([-1 1.5]);ylim([-2 2]);




%%

figure(101)
subplot(2,2,1)
hold on
plot(real(log(P1(:,i))/ST),imag(log(P1(:,i))/ST),'o')
plot([0 0],[-10 10],'LineStyle','--');
hold off
grid on
% legend('1st ','2nd ','3rd','4th','Location','northeast')
legend('Poles','Location','northeast')
title('Poles for Prop1')
xlabel('Real');ylabel('Imaginary');
xlim([-15 5]);ylim([-10 10]);
%%
subplot(2,2,2)
hold on
plot(real(log(P2(:,i))/ST),imag(log(P2(:,i))/ST),'o')
plot([0 0],[-10 10],'LineStyle','--');
hold off
grid on
% legend('1st ','2nd ','3rd','4th','Location','northeast')
legend('Poles','Location','northeast')
title('Poles for Conv1')
xlabel('Real');ylabel('Imaginary');
xlim([-15 5]);ylim([-10 10]);
%%
subplot(2,2,3)
hold on
plot(real(log(P3(:,i))/ST),imag(log(P3(:,i))/ST),'o')
plot([0 0],[-10 10],'LineStyle','--');
hold off
grid on
% legend('1st ','2nd ','3rd','4th','Location','northeast')
legend('Poles','Location','northeast')
title('Poles for Conv2')
xlabel('Real');ylabel('Imaginary');
xlim([-15 5]);ylim([-10 10]);
%% 
subplot(2,2,4)
hold on
plot(real(log(P4(:,i))/ST),imag(log(P4(:,i))/ST),'o')
plot([0 0],[-10 10],'LineStyle','--');
hold off
grid on
% legend('1st ','2nd ','3rd','4th','Location','northeast')
legend('Poles','Location','northeast')
title('Poles for Conv3')
xlabel('Real');ylabel('Imaginary');
xlim([-15 5]);ylim([-10 10]);

% subplot(2,2,[3,4])
% plot(Etime_,pose,'gx',time(1:i),BF./disps(1:i) .* mask(1:i),'ro--',time(1:i),1/lam ./ scales(1:i),'b+-.',time(1:i),X(2,1:i),'kx-')
% title('Pose Comparison')
% xlim([0 20])
% legend('Ground Truth','From Stereo Disparity','From Scaling','EKF','Location','southeast')
% grid on;
% xlabel('time [s]')
% ylabel('Estimated Velocity[m]')
% set(102,'Position',[-1500 0 800 800])

% Put current figure into Frame
% Frame1(i) = getframe(102);
%Frame2(i) = getframe(103);

% end