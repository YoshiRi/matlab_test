%% main
close all 
clear all

%% make data
% make line motion
% MakeData
% make sin motion
MakeSinData;

%% EKF estimation and show
check_cinterval = 0; % check confidence interval
rename='';
savefig = 0
% Only2D_Zupdate
% Only2D
% EKFfusion
% EKFfusion_compare
% EKFfusion_PolePlace
%  EKFfusion_jark
EKFfusion_CompareNonLinear_forpaper

%  run('./OBSapproach/OBSapproach_test.m')
%%
showComparison