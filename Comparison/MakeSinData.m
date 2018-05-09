%% config
%true depth
% sampling time 
ST = 0.033 %33ms
% ST = 0.1 %33ms
END = 20.0 % 5sec simulation

STEREO_NOISE_S = 0.5; % 1px sigma for stereo disparity noise
MONO_NOISE_S = 1; % 1?? sigma for monocular distance estimation noise

%% make ground truth data 
t = (0:ST:END).'; %time 
len = size(t,1);
freq =0.16;
phase = pi/2
% Z = 0.3 + 0.15*sin(2*pi*freq*t+pi/2);% 300mm +- 100mm
Z = 0.315 + 0.15*sin(2*pi*freq*t+phase);% 300mm +- 100mm
VZ = 2*pi*freq*0.15*cos(2*pi*freq*t+phase);
AZ = -2*pi*2*pi*freq*0.15*freq*sin(2*pi*freq*t+phase);

%% Noisy Observation 
BF = 0.065*400; % base line * focal length
StereoNoise = STEREO_NOISE_S * randn(length(t),1);
Disp = BF./Z + StereoNoise;
Zlim = 0.4;

mDisp = Disp;
INFF = 1000000000;
mDisp(mDisp>BF/Zlim) = INFF;

figure(1);
plot(t,Z)
title('GroundTruthDepth')
xlabel('time [s]')
ylabel('Object Depth [m]')
grid on;

%% Noisy Monocluar Obserbation
Z0 = 0.315 % the real distance template is taken
IMG_SIZE = 600/2% 300 times 300 pix template
MonoNoise = MONO_NOISE_S*randn(length(t),1); % sigma = 1px image noise 
Snoise = 2./(Z0./Z * IMG_SIZE).* MonoNoise;
Scale = Z0./Z + Snoise;
dScale = [1; Scale(2:length(t))./Scale(1:length(t)-1)];

hfig=figure(2)
plt = plot(t,Z,'b-.',t,BF./mDisp,'r',t,Z0./Scale,'m--')
% setfigcolor(plt,'gby')
title('Estimated Depth')
xlabel('time [s]')
ylabel('Depth [m]')
xlim([0 END/2])
legend('ground truth','3D measured','2D measured')
grid on;
    pfig = pubfig(hfig);
    pfig.LegendLoc = 'best';
    pfig.Dimension = [15 11];
    expfig(['Estimated Depth'],'-pdf');


hfig=figure(3);
plot(t,Scale,'r',t,Z0./Z,'b--')
title('Estimated Scaling')
xlabel('time [s]')
ylabel('Scaling')
legend('measured','ground truth')
grid on;
    pfig = pubfig(hfig);
    pfig.LegendLoc = 'best';
    pfig.Dimension = [15 11];
    expfig(['Estimated Scaling'],'-pdf');
