% this is newer version in Feb 6 th 2018
%
% X_prop : switch in both
% X_conv1 : switch in No stereo
% X_conv2 : switch in with stereo
% X_conv3 : no switch
%% system matrix
A = [ 1,0,0,0;0,1,ST,0;0,0,1,ST;0,0,0,1];
B = [0; ST^3/6; ST^2/2; ST];
% B = [0; 0; ST];
Q = B*B.';
% Measurment Covariance
R1 = STEREO_NOISE_S;
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

    Xmax = zeros(4,length(t));
    Xmin = zeros(4,length(t));
 

%% choose pole
p_hz = -2*pi*2;
c_pole = [p_hz,p_hz,p_hz/2 + p_hz/sqrt(3)*1i,p_hz/2 - p_hz/sqrt(3)*1i];
pole = exp(c_pole*ST);
% pole = [ 0.4408 + 0.2630i,0.4408 - 0.2630i,0.7590 + 0.0000i];
%% Xprop
BF_ = BF;
for i=2:length(t)
    % Estimate
    Xhat = A * X(:,i-1);
    Phat = A * P(:,:,i-1) * A.' + Q*Q1;

    % jacobian
    % Switch value
    if mDisp(i)== INFF
        R1_ = INFF*INFF; mD=0; R2_=R2/100;
        H2 = [0 Xhat(1) 0 0];
%         H2 = [Xhat(2) Xhat(1) 0 0];
    else
        R1_ = R1; mD = mDisp(i); R2_=R2;  
        H2 = [Xhat(2) 0 0 0];        
%         H2 = [Xhat(2) Xhat(1) 0 0];
    end
    H1 = [0 -mD^2/BF_ 0 0];
    
    % KF gain
    if mDisp(i)== INFF
        Kgain = Phat * H1.' / (H1*Phat*H1.'+R1_);
        Kgain2 = Phat * H2.' / (H2*Phat*H2.'+R2_);          
        Kpp = place(A(2:4,2:4)',A(2:4,2:4)'*[H2(2:4)]',pole(2:4))';
        Xnew =  Xhat + [Kgain [0;Kpp]]*[(mDisp(i) - BF_/Xhat(2));(1/Scale(i) - Xhat(1)*Xhat(2))];
        Pnew = (eye(4) - [Kgain [0;Kpp]]*[H1;H2])*Phat;
    else
        Kpp = place(A',A'*[H1;H2]',pole)';
        Kgain = Kpp(:,1);
        Kgain2 = Kpp(:,2);
        % update
        Xnew =  Xhat + [Kgain Kgain2]*[(mDisp(i) - BF_/Xhat(2));(1/Scale(i) - Xhat(1)*Xhat(2))];
        Pnew = (eye(4) - [Kgain Kgain2]*[H1;H2])*Phat;
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

X_prop = X;
Pole_prop = Poles;



%% Xconv1

BF_ = BF;
for i=2:length(t)
    % Estimate
    Xhat = A * X(:,i-1);
    Phat = A * P(:,:,i-1) * A.' + Q*Q1;
    
    % Switch value
    if mDisp(i)== INFF
        R1_ = INFF*INFF; mD=0; R2_=R2/100;
        H2 = [0 Xhat(1) 0 0];
        %H2 = [Xhat(2) Xhat(1) 0];
    else
        R1_ = R1; mD = mDisp(i); R2_=R2;  
        H2 = [Xhat(2) 0 0 0];        
        H2 = [Xhat(2) Xhat(1) 0 0];
    end

    
    H1 = [0 -mD^2/BF_ 0 0];
    
    % KF gain
    if mDisp(i)== INFF
        Kgain = Phat * H1.' / (H1*Phat*H1.'+R1_);
        Kgain2 = Phat * H2.' / (H2*Phat*H2.'+R2_);                
    else
        Kpp = place(A',A'*[H1;H2]',pole)';
        Kgain = Kpp(:,1);
        Kgain2 = Kpp(:,2);
    end

    % update
    Xnew =  Xhat + [Kgain Kgain2]*[(mDisp(i) - BF_/Xhat(2));(1/Scale(i) - Xhat(1)*Xhat(2))];
    Pnew = (eye(4) - [Kgain Kgain2]*[H1;H2])*Phat;

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

X_conv1 = X;
Pole_conv1 = Poles;

%% Xconv2
BF_ = BF;
for i=2:length(t)
    % Estimate
    Xhat = A * X(:,i-1);
    Phat = A * P(:,:,i-1) * A.' + Q*Q1;
    
    % Switch value
    if mDisp(i)== INFF
        R1_ = INFF*INFF; mD=0; R2_=R2/100;
        H2 = [0 Xhat(1) 0 0];
        H2 = [Xhat(2) Xhat(1) 0 0];
    else
        R1_ = R1; mD = mDisp(i); R2_=R2;  
        H2 = [Xhat(2) 0 0 0];        
        %H2 = [Xhat(2) Xhat(1) 0];
    end

    
    H1 = [0 -mD^2/BF_ 0 0];
    
    % KF gain
    if mDisp(i)== INFF
        Kgain = Phat * H1.' / (H1*Phat*H1.'+R1_);
         Kgain2 = Phat * H2.' / (H2*Phat*H2.'+R2_); 
%         Kpp = place(A',A'*[[0 1 0];H2]',pole)';
%         Kgain2 = Kpp(:,2);
    else
        Kpp = place(A',A'*[H1;H2]',pole)';
        Kgain = Kpp(:,1);
        Kgain2 = Kpp(:,2);
    end

    % update
    Xnew =  Xhat + [Kgain Kgain2]*[(mDisp(i) - BF_/Xhat(2));(1/Scale(i) - Xhat(1)*Xhat(2))];
    Pnew = (eye(4) - [Kgain Kgain2]*[H1;H2])*Phat;

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
    if mDisp(i)== INFF
        R1_ = INFF*INFF; mD=0; R2_=R2/100;
        H2 = [0 Xhat(1) 0 0];
        H2 = [Xhat(2) Xhat(1) 0 0];
    else
        R1_ = R1; mD = mDisp(i); R2_=R2;  
        H2 = [Xhat(2) 0 0 0];        
        H2 = [Xhat(2) Xhat(1) 0 0];
    end

    
    H1 = [0 -mD^2/BF_ 0 0];
    
    % KF gain
    if mDisp(i)== INFF
        Kgain = Phat * H1.' / (H1*Phat*H1.'+R1_);
        Kgain2 = Phat * H2.' / (H2*Phat*H2.'+R2_);                
    else
        Kpp = place(A',A'*[H1;H2]',pole)';
        Kgain = Kpp(:,1);
        Kgain2 = Kpp(:,2);
    end

    % update
    Xnew =  Xhat + [Kgain Kgain2]*[(mDisp(i) - BF_/Xhat(2));(1/Scale(i) - Xhat(1)*Xhat(2))];
    Pnew = (eye(4) - [Kgain Kgain2]*[H1;H2])*Phat;

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
rename = 'EKFPolePlace4'
% rename = ['EKF_1change']
% % rename = ['EKF_prop']
% showResult
%%
figure
plt=plot(t,AZ,t,X_prop(4,:),t,X_conv1(4,:),'-.',t,X_conv2(4,:),'--',t,X_conv3(4,:),':')
setfigcolor(plt,'grobp')
legend('GT','Prop','Conv1','Conv2','Conv3','Location','best')