%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Minimal working example :
%
%       -   Creates an estimator instance for a linear (runExample = 1) 
%           or non-linear (runExample = 2) kernel.
%       -   Run the kMRCD algorithm with alpha = 0.75
%           solution = kmrcd.runAlgorithm(x, alpha);
%
%   Last modified by Iwein Vranckx, 29/07/2020, 
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

    clc;
    clear all;
    close all;    
    rng(5);
    
    addpath(genpath(fileparts(which(mfilename))));
    
    color_BLUE      = [0 0.6980 0.9333];
    color_ORANGE    = [0.9333 0.4627 0];
    color_GREEN     = [0 0.6980 0.1];    
    color_GREY      = [0.25 0.25 0.25];    
    color_RED       = [0.80 0 0];
    
    %%%%    Set the example to run    
    runExample = 2;
    
    %%%%    Set the contamination degree
    epsilon = 0.2;
    
    %%%%    Set the expected amaount of regular observations.    
    alpha = 0.75;
    
    %%%%    Marker and font size
    fontSize = 10;
    mSize = 8;
    
    if runExample==2    
         N = 1000; N1 = ceil((1-epsilon)*N); N2 = N - N1;        
         data = halfkernel(N1,N2, -20, 20, 40, 5, 0.6);  
         mask = logical(data(:,3));
         x = data(:, 1:2);
         y = data(:, 3);   
         ind = randperm(size(x, 1), size(x, 1));
         x=x(ind, :);
         y=y(ind, :);
         x = rZscores(x);
         kModel = AutoRbfKernel(x); 
    else
        ndm = NewDataModel(ALYZCorrelationType(), ClusterContamination());
        [x, ~, ~,idxOutliers] = ndm.generateDataset(1000, 2, epsilon, 20);        
        y = ones(1000,1);
        y(idxOutliers) = 0;
        y = logical(y);
        x = rZscores(x);
        kModel = LinKernel(); 
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%
    %%%     Run the kMRCD algorithm... 
    %%%        
    poc = kMRCD(kModel); 
    solution = poc.runAlgorithm(x, alpha);  
    
    %   Returns:
    %   -   hsubsetIndices: indices of the h-subset elements.
    %   -   flaggedOutlierIndices:  indices of the h-flagged outliers.
    %   -   rd and cutoff    
    
    disp('We have the following solution:');
    disp(solution);
    
    rho = solution.rho;
    scfac = solution.scfac;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%
    %%%%    Visualisation    

    [rr, cc] = meshgrid(-5:0.1:5, -5:0.1:5);        
    yy=[rr(:), cc(:)];    

    Kx = kModel.compute(x(solution.hsubsetIndices, :), x(solution.hsubsetIndices, :));
    nx = size(Kx,1);
    Kt = kModel.compute(yy, x(solution.hsubsetIndices, :)); 
    Kc = Utils.center(Kx);
    Kt_c = Utils.center(Kx,Kt);
    Ktt_diag = diag(kModel.compute(yy,yy)); % Precompute
    Kxx = Ktt_diag - (2/nx)*sum(Kt,2) + (1/nx^2)*sum(sum(Kx));
    smdMesh = (1/rho)*(Kxx - (1-rho)*scfac*sum((Kt_c/((1-rho)*scfac*Kc + nx*rho*eye(nx)).*Kt_c),2)); 

    ss = logical(y);

    fig = figure(1);         
    contour(rr, cc, reshape(log(smdMesh), size(rr)), 20); hold on;
    plot(x(y>0, 1), x(y>0, 2), '.', 'color', color_GREY, 'MarkerSize', mSize); hold on;
    plot(x(y==0, 1), x(y==0, 2), '.', 'color', color_RED, 'MarkerSize', mSize);
    set(gca,'FontSize',fontSize);    
    colormap bone;
    set(gcf,'color','w');
    ylim([-4, 4]);
    title('Input dataset with marked outliers');
    saveas(fig,'images/inputdataset.png','png');
    

    fig = figure(2);         
    contour(rr, cc, reshape(log(smdMesh), size(rr)), 20); hold on;        
    plot(x(:, 1), x(:, 2), '.', 'color', color_GREY, 'MarkerSize', mSize);
    plot(x(solution.hsubsetIndices, 1), x(solution.hsubsetIndices, 2), '.', 'color', color_GREEN, 'MarkerSize', mSize);
    set(gca,'FontSize',fontSize);
    hold off;
    colormap bone;
    set(gcf,'color','w');
    ylim([-4, 4]);
    title('the h-subset');
    saveas(fig,'images/hsubset.png','png');


    fig = figure(3);        
    contour(rr, cc, reshape(log(smdMesh), size(rr)), 20); hold on;
    plot(x(:, 1), x(:, 2), '.', 'color', color_BLUE, 'MarkerSize', mSize);
    plot(x(solution.flaggedOutlierIndices, 1), x(solution.flaggedOutlierIndices, 2), '.', 'color', color_ORANGE, 'MarkerSize', mSize);
    hold off;
    colormap bone;
    set(gcf,'color','w');
    ylim([-4, 4]);
    set(gca,'FontSize',fontSize);
    title('Flagged outliers');
    saveas(fig,'images/result.png','png');

   
    
    
    
    
    
    
    
