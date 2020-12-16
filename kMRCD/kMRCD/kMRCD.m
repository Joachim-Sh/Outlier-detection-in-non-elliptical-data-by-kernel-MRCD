classdef kMRCD < handle
    %
    %kMRCD: Outlier detection in non-elliptical data by kernel MRCD, J.
    %       Schreurs, I. Vranckx et al.
    %
    %   The minimum regularized covariance determinant (MRCD) is a robust
    %   estimator for multivariate location and scatter, which detects outliers by fitting a
    %   robust covariance matrix to the data. The MRCD assumes that the observations
    %   are elliptically distributed. However, this property does not always apply to modern datasets. 
    %   Together with the time criticality of industrial processing, small n,
    %   large p problems pose a challenging problem for any data analytics procedure. Both
    %   shortcomings are solved with the proposed kernel Minimum Regularized Covariance
    %   Determinant estimator, where we exploit the kernel trick to speed-up computations.
    %   More specifically, the MRCD location and scatter matrix estimate are computed in
    %   a kernel induced feature space, where regularization ensures that the covariance matrix is well-conditioned, 
    %   invertible and defined for any dataset dimension.
    %
    %   Required input argument:
    %   x : a vector or matrix whose columns represent variables, and rows represent observations.
    %        Missing and infinite values are not allowed and should be excluded from the computations by the user.
    %   alpha : (1-alpha) measures the fraction of outliers the algorithm should
    %          resist. Any value between [0.5 <= alpha <= 1] may be specified. (default = 0.75)
    %   kModel : the kernel transformation used.
    %   
    %   Output :  The output structure 'solution' contains the final results, with the following fields :
    %       -   solution.outlyingnessIndices: outlyingness weight of each
    %       observation according.
    %       -   solution.name: the name of the best initial estimator which was used to construct the final solution.
    %       -   solution.hsubsetIndices: the h-subset element indices after
    %       C-step convergence.
    %       -   solution.obj: the MRCD objective value.
    %       -   solution.smd: the Squared Mahanobis Distance values of each
    %       observation
    %       -   solution.rho: regularisation factor used
    %       -   solution.scfac: finite sample correction factor used
    %       -   solution.rd: Mahanobis distance value of each observation.
    %       -   solution.ld: rescaled Mahanobis distance values, defined as log(0.1 + solution.rd);
    %       -   solution.cutoff: outlier flagging cut-off
    %       -   solution.flaggedOutlierIndices: indices of the flagged
    %       outliers
    %
    %   Minimal working example :
    %       -   Create an estimator instance for a linear kernel: 
    %           kmrcd = kMRCD(LinKernel());    
    %       -   Run the kMRCD algorithm with alpha = 0.75
    %           solution = kmrcd.runAlgorithm(x, alpha);
    %
    %   Last modified by Iwein Vranckx, 6/07/2020, 
    %   https://be.linkedin.com/in/ivranckx
    %   Git repository: https://github.com/ivranckx/kMRCD.git
    %   Licenced under the Non-Profit Open Software License version 3.0 (NPOSL-3.0) 
    %
    %   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
    %   INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR 
    %   PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE 
    %   FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, 
    %   TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR 
    %   THE USE OR OTHER DEALINGS IN THE SOFTWARE.
    
    
    properties (Access = private)        
        kModel;                         %   Used kernel model        
        cStepIterationsAllowed = 100;   %   Maximum number of CStep iterations allowed        
        maxcond = 50;                   %   Condition number one wants to achieve
    end
       
    methods (Access = public)
        
        function this = kMRCD(kModel)
            if ~isempty(kModel)
                this.kModel = kModel;
            else
                this.kModel = LinKernel();
            end
        end
 
         function solution = runAlgorithm(this, x, alpha)                        
            assert(alpha<=1 && alpha>=0.5, 'The percentage of regular observations, alpha, should be in [0.5-1]'); 

            %%%%    TODO:   RoS-LSSVM hack; Matlab is not type safe :/  
            %%%%            Do NOT remove the following lines - Iwein.    
            if isequal(class(this.kModel),'RbfKernel')
                disp('Warning: kMRCD switches to AutoRbfKernel in case a RbfKernel was specified!!');
                this.kModel = AutoRbfKernel(x);
            end

            K = this.kModel.compute(x, x);   
            [n, p] = size(x);
            
            %   Grab observation ranking from initial estimators
            solution = struct();            
            solution(1).outlyingnessIndices = Utils.SDO(K,alpha);
            solution(1).name = 'SDO';                
            solution(2).outlyingnessIndices = Utils.SpatialRank(K,alpha);  
            solution(2).name = 'SpatialRank';              
            solution(3).outlyingnessIndices = Utils.SpatialMedianEstimator(K,alpha);
            solution(3).name = 'SpatialMedian';
            solution(4).outlyingnessIndices = Utils.SSCM(K);
            solution(4).name = 'SSCM';
    
            scfac = Utils.MCDcons(p, alpha) ;
            
            %   For all initial estimators, do:
            rhoL = NaN(numel(solution),1);
            for index = 1:numel(solution)                
                % Compute solution for each estimator; 
                solution(index).hsubsetIndices = solution(index).outlyingnessIndices(1:ceil(n*alpha));                
                % Determine rho for each estimator            
                s = svd(Utils.center(this.kModel.compute(x(solution(index).hsubsetIndices, :), x(solution(index).hsubsetIndices, :))));
                nx = length(solution(index).hsubsetIndices);                 
                e_min = min(s); 
                e_max = max(s);
                fncond = @(rho) (nx*rho + (1-rho)*scfac*e_max)/(nx*rho + (1-rho)*scfac*e_min) - this.maxcond;
                try
                    rhoL(index) = fzero(fncond, [10^(-6),0.99]);
                catch 
                    % Find value closest to maxcond instead
                    %disp('Running grid search');
                    grid = linspace(0.000001,1-0.000001,1000);
                    objgrid = abs(arrayfun(fncond,grid));
                    rhoL(index) = min(grid(objgrid == min(objgrid)));
                end                
            end 
            
            %   Set rho as max of the rho_i's obtained for each subset in previous step                        
            rho = max(rhoL(rhoL <= max([0.1, median(rhoL)])));
            
            %   Refine each initial estimation with C-steps
            Ktt_diag = diag(K); 
            for index=1:numel(solution)                              
               for iteration = 1:this.cStepIterationsAllowed        
                    hSubset = solution(index).hsubsetIndices;    
                    Kx = this.kModel.compute(x(hSubset, :), x(hSubset, :));
                    nx = size(Kx,1);
                    Kt = this.kModel.compute(x, x(hSubset, :)); 
                    Kc = Utils.center(Kx);
                    Kt_c = Utils.center(Kx,Kt);
                    Kxx = Ktt_diag - (2/nx)*sum(Kt,2) + (1/nx^2)*sum(sum(Kx));
                    smd = (1/rho)*(Kxx - (1-rho)*scfac*sum((Kt_c/((1-rho)*scfac*Kc + nx*rho*eye(nx)).*Kt_c),2));                     
                    [~, indices] = sort(smd);                                       
                    
                    % Redefine the h-subset
                    solution(index).hsubsetIndices = indices(1:nx);                       
                    if isempty(setdiff(hSubset, solution(index).hsubsetIndices))
                        disp(['Convergence at iteration ' mat2str(iteration) ', ' solution(index).name]);                                                
                        sigma = svd(Kc);
                        sigma = (1-rho)*scfac*sigma + numel(solution(index).hsubsetIndices)*rho;                        
                        solution(index).obj =  sum(log(sigma));
                        solution(index).smd = smd;
                        break;                    
                    end                                        
               end
               assert(iteration<this.cStepIterationsAllowed, 'no C-step convergence');
            end

            %   Select the solution with the lowest objective function            
            [~, mIndex] = min([solution.obj]);     
            %   ...and remove the other solutions
            solution = solution(mIndex);                             
            disp(['-> Best estimator is ' solution.name]);
            
            solution.rho = rho;
            solution.scfac = scfac;            
            
            % Outlier flagging procedure 
            solution.rd = max(sqrt(solution.smd),0);
            solution.ld = log(0.1 + solution.rd);                                    
            [tmcd,smcd] = unimcd(solution.ld, numel(solution.hsubsetIndices));
            solution.cutoff = exp(tmcd + norminv(0.995) * smcd) - 0.1;
            solution.flaggedOutlierIndices = find(solution.rd > solution.cutoff);         
        end
        
    end
    
end

