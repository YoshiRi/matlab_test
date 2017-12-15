
%% system matrix
A = [ 1,0,0;0,1,ST;0,0,1];
B = [0; ST^2/2; ST];
% B = [0; 0; ST];
Q = B*B.';
% Measurment Covariance
R1 = STEREO_NOISE_S;
R2 = 0.0001;
Q1 = 0.1;
%% init
Pinit = diag([10000,10000,10000]);
lam0 = 1/Scale(1)/Z(1);
Xinit = [lam0 Z(1) VZ(1)];
% Xinit = [lam0 Z(1) -1];


P = zeros(3,3,length(t));
X = zeros(3,length(t));
KG = zeros(3,2,length(t));
Kpoles = zeros(3,2,length(t));
P(:,:,1) = Pinit;
X(:,1) = Xinit;

    Xmax = zeros(3,length(t));
    Xmin = zeros(3,length(t));

    
    pole =2*pi* [4,4,0] + [0,0,1];
%%
BF_ = BF;
for i=2:length(t)
    % Estimate
    Xhat = A * X(:,i-1);
    Phat = A * P(:,:,i-1) * A.' + Q*Q1;
    
    % Switch value
    if mDisp(i)== INFF
        R1_ = INFF*INFF; mD=0; 
    else
        R1_ = R1; mD = mDisp(i);  
    end
    H1 = [0 -BF_/Xhat(2)/Xhat(2) 0];
    %H1 = [0 -mD/Xhat(2) 0];
    %H1 = [0 -mD^2/BF_ 0];

    % KF gain = OBS gain
    Kgain = Phat * H1.' / (H1*Phat*H1.'+R1_);
%     Kgain  = place(A' ,H1',pole)';
    % update
    Xhat2 =  Xhat + Kgain*(mDisp(i) - BF_/Xhat(2));
    Phat2 = (eye(3) - Kgain*H1)*Phat ;

    % KF gain2
    if mDisp(i)== INFF
        H2 = [0 Xhat2(1) 0];
        R2_=R2;
    else
        H2 = [Xhat2(2) Xhat2(1) 0];  
        R2_=R2;
    end
    
    
    Kgain2 = Phat2 * H2.' / (H2*Phat2*H2.'+R2_);
%         Kgain2  = place(A' ,H2',pole)';
    % update 2
    Xnew =  Xhat2 + Kgain2*(1/Scale(i) - Xhat2(1)*Xhat2(2));
    Pnew = (eye(3) - Kgain2*H2)*Phat2;
    
    X(:,i) = Xnew;
    P(:,:,i) = Pnew;
    KG(:,1,i)=Kgain;
    KG(:,2,i)=Kgain2;
    Kpoles(:,1,i)=svd(A-Kgain*H1);
    Kpoles(:,2,i)=svd(A-Kgain2*H2);

    if check_cinterval
        [Xmax_,Xmin_]=ConfidenceInterval(Xnew,Pnew);
        Xmax(:,i) = Xmax_;
        Xmin(:,i) = Xmin_;    
    end
   
end


%%
showResultOnly

%%
figure(20)
plt = plot(t,squeeze(Kpoles(:,1,:)))


figure(21)
plt = plot(t,squeeze(Kpoles(:,2,:)),'--')
