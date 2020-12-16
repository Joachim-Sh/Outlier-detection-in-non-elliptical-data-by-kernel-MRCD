classdef AutoLaplacianKernel < handle
    
    properties (Access = public)
        sigma;
    end
    
    methods (Access = public)
        
        function this = LaplacianKernel(x)
            distances = pdist(x).^2;
            this.sigma = sqrt(median(distances));            
            disp(['AutoLaplacianKernel: Sigma = ' mat2str(this.sigma)]);            
        end
        
        function K = compute(this, Xtrain, Xtest)  
            if nargin<3
                Xtest = Xtrain;
            end            
            n=size(Xtrain, 1);    
            m=size(Xtest, 1);    
            Ka = repmat(sum(abs(Xtrain),2), 1, m);
            Kb = repmat(sum(abs(Xtest),2), 1, n);
            K = (Ka + Kb');
            K = exp(-K ./ (this.sigma));
        end
    end
    
end
