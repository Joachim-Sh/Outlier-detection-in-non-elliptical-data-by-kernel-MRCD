classdef Utils
    methods (Static)
                
        function output = MCDcons(p, alpha)
            qalpha = chi2inv(alpha, p);
            caI = gamcdf(qalpha/2, p/2 + 1,1) / alpha ;
            output = 1/caI;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% 
        %%%     Initial estimators, returns the h-subset indices in ascending
        %%%     order.
        %%% 
        
       function [gamma] = SpatialMedian(K)            
            assert(size(K, 1)==size(K, 2));
            n = size(K, 1);
            assert(n>0);
            gamma = ones(n,1)./n;
            nn = ones(n,1);
            for i = 1:10
                w = nn./sqrt(diag(K) - 2*K' * gamma + gamma'*K*gamma);
                gamma = w./sum(w);
            end                    
       end
        
       function [scale]= W_scale(x)
            
            [n,p]=size(x);
            Wc=@(x)((1-(x./4.5).^2).^2.*(abs(x)<4.5));
            sigma0=mad(x,1);
            w=Wc(x-repmat(median(x),n,1))./repmat(sigma0,n,1);
            loc=diag(x'*w)'./sum(w);
            
            rc=@(x)(min(x.^2,3^2));
            sigma0=mad(x,1);
            b=3*norminv(3/4);
            nes=n*(2*((1-b^2)*normcdf(b)-b*normpdf(b)+b^2)-1);
            scale=sigma0.^2./nes.*sum(rc((x-repmat(loc,n,1))./repmat(sigma0,n,1)));
            scale=sqrt(scale);
       end
       
       function [hIndices] = kernel_OGK(hIndices,K)                       
            n = size(K,1); 
            n_h = length(hIndices);
            K_h = K(hIndices,hIndices);
            Kt = K(:,hIndices);
           
            % Calculate Covariance matrix
            K_tilde = Utils.center(K_h);
            [U,S_F] = svd(K_tilde);
            S_F = diag(S_F);
            mask = S_F > 1000*eps; % Numerical stability
            U = U(:,mask);
            S_F = S_F(mask);
            U_scaled = U ./ repmat(sqrt(S_F'), size(U, 1), 1);

            % Step 1: Compute E and B
            o = ones(n_h, 1); 
            gamma = o./n_h;
            K_Phi_PhiTilde = (Kt -  Kt * gamma * o');
            B_F = K_Phi_PhiTilde * U_scaled; 
            lambda_F = Utils.W_scale(B_F); % Now fixed to W scale

            % Step 2: Estimate the center
            K_Adapted = K_Phi_PhiTilde*U_scaled*diag(lambda_F.^(-1))*U_scaled'*K_Phi_PhiTilde';
            gamma_c = Utils.SpatialMedian(K_Adapted);

            % Step 3: Calculate Mahalanobis
            on = ones(n, 1); 
            Kt_cCov = Kt - on * gamma_c' * Kt - Kt * gamma * o' + gamma_c' * Kt * gamma * (on * o');
            %mahal_F = diag(Kt_cCov*U_scaled*diag(lambda_F.^(-2))*U_scaled'*Kt_cCov'); % No N => scaling is in U
            mahal_F = sum((Kt_cCov*U_scaled*diag(lambda_F.^(-2))).*(Kt_cCov*U_scaled),2);
            [~,hIndices] = sort(mahal_F);
        end
        
       function [hIndices, dist, gamma] = SpatialMedianEstimator(K,h)            
            assert(size(K, 1)==size(K, 2));
            n = size(K, 1);
            
            gamma = Utils.SpatialMedian(K); 
            
            % Calculate Euclidian distance to spatial median
            dist = diag(K) - 2*sum(K.*repmat(gamma',n,1),2) + gamma'*K*gamma;
            %assert(sum((samples - gamma'*samples).^2,2)<=eps);                        
            [~, hIndices] = sort(dist); 
            
            % Kernel OGK            
            hIndices = Utils.kernel_OGK(hIndices(1:ceil(n*h)),K);
            
       end
        
        function [hIndices, gamma] = SSCM(K)
            
            assert(size(K, 1)==size(K, 2));
            % Center kernel with spatial median
            n = size(K,1);
            gamma = Utils.SpatialMedian(K); 

            % Spatial Sign Covariance matrix
            o = ones(n, 1);            
            Kc = K - o * gamma' * K - K * gamma * o' + gamma' * K * gamma * (o * o');
            D = diag((diag(K) - 2*sum(K.*repmat(gamma',n,1),2) + gamma'*K*gamma).^(-1));            
            K_tilde = D.^(1/2)*Kc*D.^(1/2);
            [U,S_F] = svd(K_tilde);
            S_F = diag(S_F);
            mask = S_F > 1000*eps; % Numerical stability
            U = U(:,mask);
            S_F = S_F(mask);
            U_scaled = U ./ repmat(sqrt(S_F'), size(U, 1), 1);
            
            %%% Kernel OGK

            % Step 1: Compute E and B
            K_Phi_PhiTilde = (K -  K * gamma * o')*D.^(1/2); % Mix Phi with Phi_tilde
            B_F = K_Phi_PhiTilde * U_scaled;  % Match
            lambda_F= Utils.W_scale(B_F); % Now fixed to W scale

            % Step 2: Estimate the center
            K_Adapted = K_Phi_PhiTilde*U_scaled*diag(lambda_F.^(-1))*U_scaled'*K_Phi_PhiTilde';
            gamma_c = Utils.SpatialMedian(K_Adapted);

            % Step 3: Calculate Mahalanobis
            K_cCov = K - o * gamma_c' * K - K * gamma * o' + gamma' * K * gamma_c * (o * o');
            K_cCov = K_cCov*D.^(1/2);
            %mahal_F = diag(K_cCov*U_scaled*diag(lambda_F.^(-2))*U_scaled'*K_cCov'); % No N => scaling is in U
            mahal_F = sum((K_cCov*U_scaled*diag(lambda_F.^(-2))).*(K_cCov*U_scaled),2);
            [~,hIndices] = sort(mahal_F);
            
        end
        
        function [hIndices, gamma] = SDO(K, h)            
            assert(size(K, 1)==size(K, 2));
            n = size(K, 1);
            gamma = zeros(n, 1);
            for direction=1:500                
                rindices = zeros(2, 1);
                while rindices(1)==rindices(2)
                    rindices  = randperm(n, 2);
                end
                lambda = zeros(n, 1);
                lambda(rindices(1))=1;
                lambda(rindices(2))=-1;
                a = K * lambda ./ sqrt( lambda' * K * lambda);                
                sdo = abs(a - median(a)) ./ mad(a);
                mask = sdo>gamma;
                gamma(mask) = sdo(mask);
            end
            
            [~, hIndices] = sort(gamma);            
            % Kernel OGK            
            hIndices = Utils.kernel_OGK(hIndices(1:ceil(n*h)),K);
        end
            
        
        function [hIndices, ook] = SpatialRank(K,h)                
            n = size(K, 1);
            ook = zeros(n,1);
            for k = 1:n                
                tmpA = K(k, k) - repmat(K(:, k), 1, n) - repmat(K(k, :), n, 1) + K;
                tmpB = sqrt(K(k,k) + diag(K)' - 2*K(k, :));
                tmpC = repmat(tmpB, n, 1) .* repmat(tmpB', 1, n);                
                mask = true(n, n);
                mask(k, :) = false;
                mask(:, k) = false;       
                ook(k) = sum(tmpA(mask) ./ tmpC(mask));                
            end
            ook = (1/n)*sqrt(ook); 
            [~, hIndices] = sort(ook);
            
            % Kernel OGK            
            hIndices = Utils.kernel_OGK(hIndices(1:ceil(n*h)),K);
        end
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% 
        %%%     Utility functions
        %%% 
        
        
        function [tmcd,smcd] = reweightedMean(y, mask)
            ncas=length(y);
            xorig=y;
            h = length(mask);
            
            initmean= mean(xorig(mask)); %initial mean
            initcov= var(xorig(mask)); %initial variance
                       
            res=(xorig-initmean).^2/initcov;
            sortres=sort(res);
            factor=sortres(h)/chi2inv(h/ncas,1);
            initcov=factor*initcov;
            res=(xorig-initmean).^2/initcov; %raw_robdist^2
            quantile=chi2inv(0.975,1);
            weights=res<=quantile; %raw-weights
            rawrd=sqrt(res);
            %reweighting procedure
            if size(weights,1)~=size(y,1)
                weights=weights';
            end
            tmcd=sum(xorig.*weights)/sum(weights);
            smcd=sqrt(sum((xorig-tmcd).^2.*weights)/(sum(weights)-1));
        end
        
        %   Centering of the kernel matrix
        function [Kc] = center(omega, Kt)
            nb_data = size(omega,1);
            Meanvec = mean(omega,2);
            MM = mean(Meanvec);
            if nargin<2
                Kc= omega-Meanvec*ones(1,nb_data)-ones(nb_data,1)*Meanvec'+MM;
            else
                nt =size(Kt,1);
                MeanvecT= mean(Kt,2);
                Kc= Kt- ones(nt,1)*Meanvec' - MeanvecT*ones(1,nb_data) +MM;
            end
        end
        
        
    end
end

