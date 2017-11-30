% Author: Yoshi Ri @ The Univ of Tokyo 
% Date : Nov 2017
% Get ConfidenceInterval of 

function [Xmax ,Xmin]= ConfidenceInterval(Xhat,Phat,p)
% Xhat : Estimated State, Phat : Estimated Covariance, p : percentage of
% prediction

if nargin==2
p=0.95;
%print('set default p=95%!')
end

% Get dimension
Dim = size(Xhat,1);
alpha = p;
bound = chi2inv(alpha,Dim);

% Set symbolic equation
Pinv= inv(Phat);

% Get upper and lower boundary
Xmax=zeros(size(Xhat));
Xmin=zeros(size(Xhat));
for i=1:Dim
    Xmin(i)=Xhat(i)- sqrt(bound/Pinv(i,i));
    Xmax(i)=Xhat(i)+ sqrt(bound/Pinv(i,i));
end

end
