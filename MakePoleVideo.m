%% Repetitive


% make video
Frate = 3;FLen = 102;
Frame1(FLen) = struct('cdata',[],'colormap',[]);
Frame2(FLen) = struct('cdata',[],'colormap',[]);
figure(102)
for i = 1:FLen
close(102)
figure(102)
subplot(2,2,1)
hold on
viscircles([0 0],1,'LineStyle','--');
plot(real(P1(:,i)),imag(P1(:,i)),'o')
hold off
grid on
legend('1st ','2nd ','3rd','4th','Location','northeast')
title('Poles for stereo observation')
xlabel('Real');ylabel('Imaginary');
xlim([-1 1.5]);ylim([-2 2]);
subplot(2,2,2)
hold on
viscircles([0 0],1,'LineStyle','--');
plot(real(P2(:,i)),imag(P2(:,i)),'o')
hold off
grid on
legend('1st ','2nd ','3rd','4th','Location','northeast')
title('Poles for mono observation')
xlabel('Real');ylabel('Imaginary');
xlim([-1 1.5]);ylim([-2 2]); 
subplot(2,2,[3,4])
plot(Etime_,pose,'gx',time(1:i),BF./disps(1:i) .* mask(1:i),'ro--',time(1:i),1/lam ./ scales(1:i),'b+-.',time(1:i),X(2,1:i),'kx-')
title('Pose Comparison')
xlim([0 20])
legend('Ground Truth','From Stereo Disparity','From Scaling','EKF','Location','southeast')
grid on;
xlabel('time [s]')
ylabel('Estimated Velocity[m]')
set(102,'Position',[-1500 0 800 800])

% Put current figure into Frame
Frame1(i) = getframe(102);
%Frame2(i) = getframe(103);

end

% % convert the Frame to movie and show
% figure(2);
% movie(Frame,1);

%% write to video
v1 = VideoWriter('pole12.avi');
% v2 = VideoWriter('pole2.avi');
v1.FrameRate = Frate; % Framerate
% v2.FrameRate = Frate; % Framerate
open(v1);
writeVideo(v1,Frame1);
close(v1);
% open(v2);
% writeVideo(v2,Frame2);
% close(v2);

