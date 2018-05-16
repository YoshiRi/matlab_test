%% main
close all 
clear all

%% Simulation Config

ST = 0.033 %33ms
% ST = 0.1 %33ms
END = 20.0 % 5sec simulation

STEREO_NOISE_S = 0.5; % 1px sigma for stereo disparity noise
MONO_NOISE_S = 1; % 1?? sigma for monocular distance estimation noise

%% make data
MakeSinData;

%% EKF estimation and show
check_cinterval = 0; % check confidence interval
rename='';
% Only2D_Zupdate
% Only2D
% EKFfusion
EKFfusion_compare
% EKFfusion_PolePlace
%  EKFfusion_jark
%  run('./OBSapproach/OBSapproach_test.m')
%%
showComparison