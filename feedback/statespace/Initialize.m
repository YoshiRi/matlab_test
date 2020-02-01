Xlm = zeros(3,length(t));
Dists = zeros(2,length(t));
U = zeros(1,length(t));

Pinit = diag([10000,10000,10000]);
lam0 = 1/dref;
Xinit = [lam0*0.6 dref 0];
% Xinit = [lam0 Z(1) -1];

% proposed
P = zeros(3,3,length(t));
X = zeros(3,length(t));
KG = zeros(3,2,length(t));
Poles = zeros(3,length(t));
P(:,:,1) = Pinit;
X(:,1) = Xinit;

Xmax = zeros(3,length(t));
Xmin = zeros(3,length(t));

Xlm(:,1) = [-dref;vf0;0];