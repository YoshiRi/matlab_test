%% 
% this is newer version in Feb 6 th 2018
%
% X_prop : switch in both
% X_conv1 : switch in No stereo
% X_conv2 : switch in with stereo
% X_conv3 : no switch
%% ACC noise
ACC_NOISE_S = 0.001;
ACCNoise = ACC_NOISE_S * randn(length(t),1);
AZ_ = AZ + ACCNoise;
% Z+Lambda
% ZpL = Z+lam0+ACCNoise;
%% choose pole
p_hz = -2*pi*2;
c_pole = [p_hz,p_hz,p_hz/2 + p_hz/sqrt(3)*1i,p_hz/2 - p_hz/sqrt(3)*1i];
pole = exp(c_pole*ST);

%% system matrix
A = [ 1,0,0,0;0,1,ST,0;0,0,1,ST;0,0,0,1];
B = [0; ST^3/6; ST^2/2; ST];
% B = [0; 0; ST];
Q = B*B.';
% Measurment Covariance
R1 = ACC_NOISE_S;
R2 = 0.01;
Q1 = 0.1;
%% init
Pinit = diag([10000,10000,10000,10000]);
lam0 = 1/Scale(1)/Z(1);
Xinit = [lam0*0.6 Z(1) VZ(1) AZ(1)];
% Xinit = [lam0 Z(1) -1];

% proposed
P = zeros(4,4,length(t));
X = zeros(4,length(t));
KG = zeros(4,2,length(t));
Poles = zeros(4,length(t));
P(:,:,1) = Pinit;
X(:,1) = Xinit;

Xmax = zeros(3,length(t));
Xmin = zeros(3,length(t));

H1 = [0 0 0 1];
H2 = [0 0 0 1];

pH1 = H1;
pH2 = H2;

%% Xprop
BF_ = BF;
for i=2:length(t)
    % Estimate
    Xhat = A * X(:,i-1);
    Phat = A * P(:,:,i-1) * A.' + Q*Q1;
    
    % Switch value
    R1_ = R1; mD=0; R2_=R2;
    H1 = [0 0 0 1];
    if mod(i,12)==0
        H2 = [Xhat(2) Xhat(1) 0 0];
    % KF gain
    Kgain = Phat * H1.' / (H1*Phat*H1.'+R1_);
    Kgain2 = Phat * H2.' / (H2*Phat*H2.'+R2_);

    % update
    Xnew =  Xhat + [Kgain Kgain2]*[(AZ_(i) - Xhat(4));(1/Scale(i) - Xhat(1)*Xhat(2))];
    Pnew = (eye(4) - [Kgain Kgain2]*[H1;H2])*Phat;

    else
    % KF gain
    Kgain = Phat * H1.' / (H1*Phat*H1.'+R1_);
    Kgain2 = Kgain;

    % update
    Xnew =  Xhat + [Kgain]*[(AZ_(i) - Xhat(4))];
    Pnew = (eye(4) - [Kgain ]*[H1])*Phat;
    end
    
    % update 2
    X(:,i) = Xnew;
    P(:,:,i) = Pnew;
    KG(:,1,i)=Kgain;
    KG(:,2,i)=Kgain2;
    Poles(:,i) = eig((A-[Kgain Kgain2]*[H1;H2]*A)*(A-[KG(:,1,i-1) KG(:,2,i-1)]*[pH1;pH2]*A));
    pH1 = H1;
    pH2 = H2;
    
    if check_cinterval
        [Xmax_,Xmin_]=ConfidenceInterval(Xnew,Pnew);
        Xmax(:,i) = Xmax_;
        Xmin(:,i) = Xmin_;    
    end
   
end

X_prop = X;
Pole_prop = Poles;
%% Xconv1
BF_ = BF;
for i=2:length(t)
    % Estimate
    Xhat = A * X(:,i-1);
    Phat = A * P(:,:,i-1) * A.' + Q*Q1;
    
    % Switch value
    R1_ = R1; mD=0; R2_=R2;
    H1 = [0 0 0 1];
    if mod(i,4)==0
        H2 = [Xhat(2) Xhat(1) 0 0];
    % KF gain
    Kgain = Phat * H1.' / (H1*Phat*H1.'+R1_);
    Kgain2 = Phat * H2.' / (H2*Phat*H2.'+R2_);

    % update
    Xnew =  Xhat + [Kgain Kgain2]*[(AZ_(i) - Xhat(4));(1/Scale(i) - Xhat(1)*Xhat(2))];
    Pnew = (eye(4) - [Kgain Kgain2]*[H1;H2])*Phat;

    else
    % KF gain
    Kgain = Phat * H1.' / (H1*Phat*H1.'+R1_);

    % update
    Xnew =  Xhat + [Kgain]*[(AZ_(i) - Xhat(4))];
    Pnew = (eye(4) - [Kgain ]*[H1])*Phat;
    end

    % update 2
    X(:,i) = Xnew;
    P(:,:,i) = Pnew;
    KG(:,1,i)=Kgain;
    KG(:,2,i)=Kgain2;
    Poles(:,i) = eig((A-[Kgain Kgain2]*[H1;H2]*A)*(A-[KG(:,1,i-1) KG(:,2,i-1)]*[pH1;pH2]*A));
    pH1 = H1;
    pH2 = H2;

    if check_cinterval
        [Xmax_,Xmin_]=ConfidenceInterval(Xnew,Pnew);
        Xmax(:,i) = Xmax_;
        Xmin(:,i) = Xmin_;    
    end
   
end

X_conv1 = X;
Pole_conv1 = Poles;

%% Xconv2
BF_ = BF;
for i=2:length(t)
    % Estimate
    Xhat = A * X(:,i-1);
    Phat = A * P(:,:,i-1) * A.' + Q*Q1;
    
    % Switch value
    R1_ = R1; mD=0; R2_=R2;
    H1 = [0 0 0 1];
    if mod(i,6)==0
        H2 = [Xhat(2) Xhat(1) 0 0];
    % KF gain
    Kgain = Phat * H1.' / (H1*Phat*H1.'+R1_);
    Kgain2 = Phat * H2.' / (H2*Phat*H2.'+R2_);

    % update
    Xnew =  Xhat + [Kgain Kgain2]*[(AZ_(i) - Xhat(4));(1/Scale(i) - Xhat(1)*Xhat(2))];
    Pnew = (eye(4) - [Kgain Kgain2]*[H1;H2])*Phat;

    else
    % KF gain
    Kgain = Phat * H1.' / (H1*Phat*H1.'+R1_);

    % update
    Xnew =  Xhat + [Kgain]*[(AZ_(i) - Xhat(4))];
    Pnew = (eye(4) - [Kgain ]*[H1])*Phat;
    end

    % update 2
    X(:,i) = Xnew;
    P(:,:,i) = Pnew;
    KG(:,1,i)=Kgain;
    KG(:,2,i)=Kgain2;
    Poles(:,i) = eig(A-[Kgain Kgain2]*[H1;H2]*A);

    if check_cinterval
        [Xmax_,Xmin_]=ConfidenceInterval(Xnew,Pnew);
        Xmax(:,i) = Xmax_;
        Xmin(:,i) = Xmin_;    
    end
   
end

X_conv2 = X;
Pole_conv2 = Poles;

%% Xconv3
BF_ = BF;
for i=2:length(t)
    % Estimate
    Xhat = A * X(:,i-1);
    Phat = A * P(:,:,i-1) * A.' + Q*Q1;
    
    % Switch value
    R1_ = R1; mD=0; R2_=R2;
    H1 = [0 0 0 1];
    if mod(i,2)==0
        H2 = [Xhat(2) Xhat(1) 0 0];
    % KF gain
    Kgain = Phat * H1.' / (H1*Phat*H1.'+R1_);
    Kgain2 = Phat * H2.' / (H2*Phat*H2.'+R2_);

    % update
    Xnew =  Xhat + [Kgain Kgain2]*[(AZ_(i) - Xhat(4));(1/Scale(i) - Xhat(1)*Xhat(2))];
    Pnew = (eye(4) - [Kgain Kgain2]*[H1;H2])*Phat;

    else
    % KF gain
    Kgain = Phat * H1.' / (H1*Phat*H1.'+R1_);

    % update
    Xnew =  Xhat + [Kgain]*[(AZ_(i) - Xhat(4))];
    Pnew = (eye(4) - [Kgain ]*[H1])*Phat;
    end

    % update 2
    X(:,i) = Xnew;
    P(:,:,i) = Pnew;
    KG(:,1,i)=Kgain;
    KG(:,2,i)=Kgain2;
    Poles(:,i) = eig(A-[Kgain Kgain2]*[H1;H2]*A);

    if check_cinterval
        [Xmax_,Xmin_]=ConfidenceInterval(Xnew,Pnew);
        Xmax(:,i) = Xmax_;
        Xmin(:,i) = Xmin_;    
    end
   
end

X_conv3 = X;
Pole_conv3 = Poles;

%%
rename = 'debug'
% rename = ['EKF_1change']
% % rename = ['EKF_prop']
% showResult

showPole

figure(30)
plt=plot(t,AZ,t,X_prop(4,:),t,X_conv1(4,:),'-.',t,X_conv2(4,:),'--',t,X_conv3(4,:),':')
setfigcolor(plt,'grobp')
legend('GT','Prop','Conv1','Conv2','Conv3','Location','best')

showComparison
