
%%
ST = 0.033 %33ms
% ST = 0.1 %33ms
END = 30.0 % 20sec simulation

STEREO_NOISE_S = 1.5; % 1px sigma for stereo disparity noise
MONO_NOISE_S = 0.1; % 1?? sigma for monocular distance estimation noise

%% sensor
BF = 0.065*400*10; % base line * focal length
Zlim = 0.3;
INFF = 1000000000;

%% make forward vehicle data 
t = (0:ST:END).'; %time 
len = size(t,1);

vf0 = 40*1000/3600; %40km/h to m/s
vf1 = 60*1000/3600;
Tst = 150;
Tfi = 300;

vfront = vf0 * ones(len,1);
vfront(Tst+1:Tfi)=(1:(Tfi-Tst))/(Tfi-Tst)*(vf1-vf0)+vf0;
vfront(Tfi+1:end)=vf1;

pfront = cumsum(vfront*ST);
figure(1)
subplot(121)
plot(t,vfront)
subplot(122)
plot(t,pfront)


%% Make noise
STEREO_NOISE_S = 0.5; % 1px sigma for stereo disparity noise
MONO_NOISE_S = 0.1; % 1?? sigma for monocular distance estimation noise

dref = 20;
Z0 = dref;

IMG_SIZE = 600/2% 300 times 300 pix template
MonoNoise = MONO_NOISE_S*randn(length(t),1); % sigma = 1px image noise 
Snoise = 2./(Z0 * IMG_SIZE).* MonoNoise;
StereoNoise = STEREO_NOISE_S * randn(length(t),1);

