classdef PolyKernel < handle
    
    properties (Access = private)
        degree;
    end
    
    methods (Access = public)
        
        function this = PolyKernel(degree)
            this.degree = degree;
        end
        
        function K = compute(this, Xtrain, Xtest)                 
            if nargin<3
                Xtest = Xtrain;
            end  
            K = (Xtrain * Xtest' + 1).^this.degree;
        end
    end
    
end

