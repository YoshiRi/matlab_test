% 

figure(9)
plot([t(1) t(end)],[1/Z0 1/Z0],'--',t,X(1,:))
grid on;
title('Mag Estimation')

figure(10)
plot(t,Z,t,Xlm(1,:),'--')
legend('obj','self');
grid on;
title('Obj and Self Position')

figure(11)
plot(t,Dists(1,:),t,X(2,:),'--')
title('Distance Est')
legend('Ground Truth','Estimation')
grid on;

figure(12)
plot(t,Dists(2,:),t,X(3,:),'--')
title('Velocity Est')
legend('Ground Truth','Estimation')
grid on;

figure(13)
plot(t,Dists(1,:),[t(1) t(end)],[Dref Dref],'--')
grid on;
legend('Tracked','Ref')
title('Distance reference and tracked result')


figure(14)
plot(t,U);
grid on;
title('Input amplitude')

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
plt2 = plot(real(eig(FB)),imag(eig(FB)),'bx');
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