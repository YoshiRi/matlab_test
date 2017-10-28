clear all
close all

set(0, 'DefaultLineLineWidth', 2) 

%% config
%true depth
% sampling time 
ST = 0.033 %33ms
END = 5.0 % 5sec simulation
STEREO_NOISE_S = 1; % 1px sigma for stereo disparity noise
MONO_NOISE_S = 1; % 1?? sigma for monocular distance estimation noise

%% make ground truth data 
t = (0:ST:END).'; %time 
len = size(t,1);
freq =0.4;
Z = 0.3 + 0.1*sin(2*pi*freq*t);% 300mm +- 100mm
VZ = 2*pi*freq*0.1*cos(2*pi*freq*t);
Z = 0.8 - 0.1*t;% 300mm +- 100mm
VZ = -0.1+0*t;

%% Noisy Observation 
BF = 0.065*400; % base line * focal length
StereoNoise = STEREO_NOISE_S * randn(length(t),1);
Disp = BF./Z + StereoNoise;


mDisp = Disp;
INFF = 1000000000;
mDisp(mDisp>BF/0.275) = INFF;

figure(1);
plot(t,Z)
title('GroundTruthDepth')
xlabel('time [s]')
ylabel('Object Depth [m]')
grid on;
figure(2)
plot(t,BF./mDisp,'r',t,Z,'b--')
title('Estimated Depth')
xlabel('time [s]')
ylabel('Depth [m]')
legend('measured','ground truth')
grid on;

%% Noisy Monocluar Obserbation
Z0 = 0.3 % the real distance template is taken
IMG_SIZE = 600/2% 300 times 300 pix template
MonoNoise = MONO_NOISE_S*randn(length(t),1); % sigma = 1px image noise 
Snoise = 2./(Z0./Z * IMG_SIZE).* MonoNoise;
Scale = Z0./Z + Snoise;

figure(3);
plot(t,Scale,'r',t,Z0./Z,'b--')
title('Estimated Scaling')
xlabel('time [s]')
ylabel('Scaling')
legend('measured','ground truth')
grid on;
