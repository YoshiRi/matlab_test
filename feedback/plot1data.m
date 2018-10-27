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
