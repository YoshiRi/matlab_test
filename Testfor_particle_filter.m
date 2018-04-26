%%  test for particle filter
% Problem Setting is here
% https://ja.wikipedia.org/wiki/%E3%82%AB%E3%83%AB%E3%83%9E%E3%83%B3%E3%83%95%E3%82%A3%E3%83%AB%E3%82%BF%E3%83%BC#%E5%BF%9C%E7%94%A8%E4%BE%8B

% settings
dt = 0.001; % time sample
s_sys = 10000; % system noise
s_out = 0.05; % measurement noise
s_rout = 0.1;
% system eq
F = [1,dt;0,1];
G = [dt^2/2;dt];
covG = s_sys^2*G*G';
% output eq
H = [1 0];
R = s_out^2;
%init val
Xini = [0 ;0];

% state_equation
stateEq = @(x,u) F*x
outEq = @(x) H*x
% Init particle filter
pfilt = particle_filter(100,Xini,stateEq,outEq,covG,R);

%% Make Reference data
time = 0:100;
time = time' * dt;
input = 3*sin(2*pi*20*time); % sin speed input

Xref = zeros(2,size(time,1));
Xest = Xref;
Yref = zeros(1,size(time,1));
for i = 2:size(time,1)
    Xref(:,i) = stateEq(Xref(:,i-1)) + [0;input(i-1)];
    Yref(:,i) = outEq(Xref(:,i));
end
% noisy obserbation
Yobs = Yref + s_rout*randn(1,size(time,1));
Vobs = diff(Yobs)/dt;
%% show reference
figure(1)
plot(time,Xref)
legend('X','V')
figure(2)
plot(time,Yref,'--',time,Yobs,'-')
legend('Xactual','Xobs')

%% 
 for i = 2:size(time,1)
     Xest(:,i) = pfilt.predict(0,Yobs(i));
 end

%% 
figure(3)
plot(time,Yref,'--',time,Yobs,'-',time,Xest(1,:))
legend('Xactual','Xobs','Xest')
figure(4)
plot(time,Xref(2,:),'--',time,Xest(2,:),time(2:end),Vobs)
legend('Vactual','Vest','Raw Derivative')

figure(5)
plot(Xref(1,end),Xref(2,end),'o',Xest(1,end),Xest(2,end),'x',pfilt.particle_state(1,:),pfilt.particle_state(2,:),'+')
legend('Xactual','Xaverage','Xparticles')