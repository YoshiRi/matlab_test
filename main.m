%% main
close all 
clear all

%% make data
% make line motion
%MakeData
% make sin motion
MakeSinData;

%% EKF estimation and show
check_cinterval = 1; % check confidence interval
rename='';
% Only2D_Zupdate
%Only2D
%  EKFfusion
 EKFfusion_jark
 