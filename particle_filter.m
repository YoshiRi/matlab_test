classdef particle_filter < handle
% Initialization:
%    obj = particle_filter(n,Xini,stateEq,outputEq,SsysNoise,SoutputNoise)
%
% Usage:
%    Xest = obj.predict(input,output)
%
% Public parameter:
%     particle_state, particle_weight

    properties (SetAccess = protected)
        % Estimated Value
        pnum = 1;
        dim = 1;
        % Noise
        sigma_sys = 0; % diag([1,2,3])
        sigma_out = 0; % diag([1,2,3])
        % Function for updating: Example
        state_equation = []; % @(x,u) A*x+B*u
        output_equation = []; % @(x,u) C*x+D*u
    end
    properties (SetAccess = public)
        particle_state = [];
        particle_pstate = [];
        particle_weight = [];
        particle_pweight = [];
    end
    
    methods
        
        % constructer
        function obj = particle_filter(n,Xini,stateEq,outputEq,SsysNoise,SoutputNoise)
           obj.pnum = n;
           obj.dim = size(Xini,1);
           obj.state_equation = stateEq;
           obj.output_equation = outputEq;
           obj.sigma_sys = SsysNoise;
           obj.sigma_out = SoutputNoise;
           obj.init_particle(Xini);
        end
        
        function init_particle(obj,Xini)
            obj.particle_state = repmat(Xini,1,obj.pnum)+sqrt(obj.sigma_sys)*randn(obj.dim,obj.pnum);
            obj.particle_pstate = obj.particle_state;
            obj.particle_weight = zeros(1,obj.pnum) + 1/obj.pnum;
            obj.particle_pweight = obj.particle_weight;
        end
        
        % Prediction and Update
        function Xest = predict(obj,input,output)
            % for loop is not desirable for faster implementation
            for i = 1:obj.pnum
                % state eq
                obj.particle_pstate(:,i) = obj.state_equation(obj.particle_state(:,i),input)+sqrt(obj.sigma_sys)*randn(obj.dim,1);
                % output eq
                dz = output - obj.output_equation(obj.particle_pstate(:,i));
                obj.particle_pweight(:,i) = obj.particle_weight(:,i) * obj.Gauss(dz,0,obj.sigma_out);
            end
            lambda = sum(obj.particle_pweight);
            obj.particle_pweight = obj.particle_pweight/lambda;
            obj.Resampling();
            Xest = obj.particle_state * obj.particle_weight';
        end
        
        
        % For multiple obserbation based filtering
%         function state_transition(obj,input,output)
%             % For only linear equation
%             Minput = repmat(input,1,obj.pnum);
%             Mnoise = sqrt(obj.sigma_sys)*randn(obj.dim,obj.pnum);
%             obj.particle_pstate = obj.state_equation(obj.particle_state,Minput)+Mnoise;
%             dz1 = output()
%             dz2 = 
%             for i = 1:obj.pnum
%                 % state eq
%                 obj.particle_pstate(:,i) = obj.state_equation(obj.particle_state(:,i),input)+sqrt(obj.sigma_sys)*randn(obj.dim,1);
%                 % output eq
%                 dz = output - obj.output_equation(obj.particle_pstate(:,i));
%                 obj.particle_pweight(:,i) = obj.particle_weight(:,i) * obj.Gauss(dz,0,obj.sigma_out);
%             end
%             lambda = sum(obj.particle_pweight);
%             obj.particle_pweight = obj.particle_pweight/lambda;
%             obj.Resampling();
%             Xest = obj.particle_state * obj.particle_weight';
%         end
       
        
        function p=Gauss(obj,x,mu,sigma2)
            p = 1/sqrt(2*pi)^obj.dim / sqrt(det(sigma2)) * exp(-0.5*(x-mu).'/sigma2*(x-mu));
        end
        
        function Resampling(obj)
            %リサンプリングを実施する関数
            pw = obj.particle_pweight;
            NTh = obj.pnum/2;
            NP = obj.pnum;

            %アルゴリズムはLow Variance Sampling
            Neff=1.0/(pw*pw');
            if Neff<NTh %リサンプリング
                wcum=cumsum(pw);
                base=cumsum(pw*0+1/NP)-1/NP;%乱数を加える前のbase
                resampleID=base+rand/NP;%ルーレットを乱数分増やす
                ppx =obj.particle_pstate;%データ格納用
                ind=1;%新しいID
                for ip=1:NP
                    while(resampleID(ip)>wcum(ind))
                        ind=ind+1;
                    end
                    obj.particle_state(:,ip)=ppx(:,ind);%LVSで選ばれたパーティクルに置き換え
                    obj.particle_weight(ip)=1/NP;%尤度は初期化
                end
            end
            
        end
    end
    
end