% 

figure(9)
plot([t(1) t(end)],[1/Z0 1/Z0],'--',t,X(1,:))
grid on;
title('Magnitude Estimation')
xlabel('time [s]')
ylabel('Magnitude constant')

figure(10)
plot(t,pfront,t,Xlm(1,:),'--')
legend('Front','Back','Location','best');
grid on;
title('Front and back car position')
xlabel('time [s]')
ylabel('position [m]')

figure(11)
plot(t,Dists(1,:),t,X(2,:),'--')
title('Distance estimation')
legend('Ground Truth','Estimation','Location','best')
grid on;
xlabel('time [s]')
ylabel('relative position [m]')

figure(12)
plot(t,vfront,t,Xlm(2,:),'--')
title('Velocity tracking results')
legend('Front','Back','Location','best')
grid on;
xlabel('time [s]')
ylabel('velocity [m]')


figure(13)
plot(t,Dists(2,:),t,X(3,:),'--')
title('Velocity estimation')
legend('Ground Truth','Estimation','Location','best')
grid on;
xlabel('time [s]')
ylabel('relative velocity [m/s]')


figure(14)
plot([t(1) t(end)],[dref dref],t,Dists(1,:),'--')
grid on;
legend('Reference','Tracked','Location','best')
title('Distance reference and tracked result')
xlabel('time [s]')
ylabel('relative position [m]')


figure(15)
plot(t,U);
grid on;
title('Inputs acceleration for back car')
xlabel('time [s]')
ylabel('acceleration [m/s^2]')
%%
Pvals = zeros(3,length(t));
for i = 1:length(t)
  Pvals(:,i) = 1./diag(inv( squeeze(P(:,:,i)) ));
end

figure(20)
plot(t,Pvals)
ylim([0 1])
legend('Mag','Pos','Vel')
%%
figure(100)
polarplot(Poles,'o');

figure(101)
for i = 1:length(t)
plt =plot(real(log(Poles(:,i))/ST),imag(log(Poles(:,i))/ST),'o');
plt.Color = [1,1-i/length(t),i/length(t)];
hold on
end
plt2 = plot(real(PoleC),imag(PoleC),'bx');
legend([plt,plt2],{'poles for EKF', 'poles for Control'})
hold off
xlim([-3 1])
grid on
% legend('poles for EKF')

%% figure

figure(200)

Covs =  zeros(length(t),3); 
for i = 1:length(t)
    Covs(i,:) = diag(squeeze(P(:,:,i)));
end

plot(t,Covs)
grid on 
legend('\sigma_\lambda','\sigma_D','\sigma_V')