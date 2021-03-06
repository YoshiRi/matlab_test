
%% system matrix
A = [ 1,0,0;0,1,ST;0,0,1];
B = [0; ST^2/2; ST];
Q = B*B.';
% Measurment Covariance
R1 = STEREO_NOISE_S;
R2 = 0.005;

%% init
Pinit = diag([100,100,100]);
lam0 = 1/Z0;
Xinit = [lam0 Z(1) VZ(1)];

P = zeros(3,3,length(t));
X = zeros(3,length(t));
KG = zeros(3,2,length(t));
P(:,:,1) = Pinit;
X(:,1) = Xinit;

%%
BF_ = BF;
for i=2:length(t)
    % Estimate
    Xhat = A * X(:,i-1);
    Phat = A*P(:,:,i-1)*A.' + Q;

    
    
    
    % Switch value
    if mDisp(i)== INFF
        R1_ = INFF*INFF; mD=0;
    else
        R1_ = R1; mD = mDisp(i);        
    end
    %H1 = [0 -BF_/Xhat(2)/Xhat(2) 0];
    %H1 = [0 -mD/Xhat(2) 0];
    H1 = [0 -mD^2/BF_ 0];

    H2 = [Xhat(2) Xhat(1) 0];
    H2 = [0 1/Z0 0];
    
    Kgain = Phat * H2.' / (H2*Phat*H2.'+R2);
    % update 2
    Xhat2 =  Xhat + Kgain*(1/Scale(i) - Xhat(1)*Xhat(2));
    Phat2 = (eye(3) - Kgain*H2)*Phat;

    
    % KF gain
    Kgain2 = Phat2 * H1.' / (H1*Phat2*H1.'+R1_);
    % update
    Xnew =  Xhat2 + Kgain2*(mDisp(i) - BF_/Xhat2(2));
    Pnew = (eye(3) - Kgain2*H1)*Phat2 ;

    % KF gain2
%     if mDisp(i)== INFF
%         H2 = [0 Xhat2(1) 0];
%     else
%         H2 = [Xhat2(2) Xhat2(1) 0];
%     end
    X(:,i) = Xnew;
    P(:,:,i) = Pnew;
    KG(:,1,i)=Kgain;
    KG(:,2,i)=Kgain2;
end


%% figure
showResult