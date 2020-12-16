classdef RbfKernel < handle
    
    properties (Access = public)
        sigma;
    end
    
    methods (Access = public)
        
        function this = RbfKernel(bandwidth)
            this.sigma = bandwidth;
        end
        
        function updateKernel(this, bandwidth)
           this.sigma = bandwidth; 
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

