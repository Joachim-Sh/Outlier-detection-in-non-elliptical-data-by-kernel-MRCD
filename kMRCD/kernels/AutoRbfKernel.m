classdef AutoRbfKernel < handle
    
    properties (Access = public)
        sigma;
    end
    
    methods (Access = public)
        
        function this = AutoRbfKernel(x)
            distances = pdist(x).^2;
            this.sigma = sqrt(median(distances));            
            disp(['AutoRbfKernel: Sigma = ' mat2str(this.sigma)]);
        end
        
        function K = compute(this, Xtrain, Xtest)  
            if nargin<3
                Xtest = Xtrain;
            end            
            n=size(Xtrain, 1);    
            m=size(Xtest, 1);    
            Ka = repmat(sum(Xtrain.^2,2), 1, m);
            Kb = repmat(sum(Xtest.^2,2), 1, n);
            K = (Ka + Kb' - 2 .* (Xtrain * Xtest'));
            K = exp(-K ./ (2* this.sigma^2));
        end
    end
    
end

