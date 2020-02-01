%% main
close all 
clear all

%% make data
% make line motion
%MakeData
% make sin motion
MakeSinData;

%% EKF estimation and show
check_cinterval = 0; % check confidence interval
rename='';

%OBSapproach_test
Poleplacement
%  EKFfusion_jark
%  run('./OBSapproach/OBSapproach_test.m')
%%
showComparison