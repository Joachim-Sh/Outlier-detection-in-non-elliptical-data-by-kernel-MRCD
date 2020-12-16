classdef ALYZCorrelationType < handle
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        
        function cor = covariance2correlationInPlace(~, covariance)	
            D = diag(diag(covariance).^(-1/2));
            cor= D*covariance*D;
        end
        
		%	Generate random correlation matrix with given condition number
        function corEstimation = generateCorr(this, p)		
			conditionNumber = 100;			
			HI = conditionNumber; % set HI and LO according to your problem.
			LO = 1.0;
			aRandomVector = sort(LO + (HI-LO).*rand(p-2,1));
			L = diag([1; aRandomVector; conditionNumber]);
            
			Y  = randn(p, p);		
            [U,~,~]=svd(Y' * Y);
			
			corEstimation = (U * L * U');
			for iteration = 1:100
                
				corEstimation = this.covariance2correlationInPlace(corEstimation);
                [U,S,~] = svd(corEstimation);
				L = S;
 				oldConditionNumber = L(1, 1) / L(p, p);
				
 				if (abs(conditionNumber - oldConditionNumber) < 10^(-10))
                    %display(['Condition number is ' mat2str(conditionNumber)]);
 					break;				
                end
				
				L(1, 1) = conditionNumber * L(p, p);
				corEstimation = (U * L * U');
            end

			corEstimation = this.covariance2correlationInPlace(corEstimation);
        end

        function mu = generateLocation(~, p)		
			mu = zeros(p, 1);
        end
    
    end
    
end

