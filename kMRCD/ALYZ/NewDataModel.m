classdef NewDataModel < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        CorrelationType;
        ContaminationType;
    end
    
    methods
        
        function this = NewDataModel(CorrelationType, ContaminationType)
            this.CorrelationType = CorrelationType;
            this.ContaminationType = ContaminationType;
        end
        
        function Sigma = generateA09cormatrix(~,d, correl) 
              % generates A09 correlation matrix; correl = 0.9 is the default choice
              % Note: -correl is used as base value (!)
              columns = repmat((1:d),d,1);
              rows    = repmat((1:d)',1,d);
              Sigma   = -correl*ones(d);
              Sigma   = Sigma.^(abs(columns - rows)) ;
        end

        
        function [samples, tCorrelation, tLocation,cindices] = generateDataset(this, n,  p, eps, k)	
            contaminationDegree = floor(n * eps);

            %	Generate correlation matrix
            tCorrelation = this.CorrelationType.generateCorr(p);
            tLocation = this.CorrelationType.generateLocation(p);

            %	Generate multivariate normal data with location zero and correlation;
            samples = mvnrnd(tLocation, tCorrelation, n);

            %   Generate contamination 
            [U,~,~] = svd(tCorrelation);
            replacement = U(:,end);
            
            delta = replacement - tLocation;
            smd = sqrt((delta' * inv(tCorrelation) * delta));
            replacement = replacement * (k/smd); 
            
            % Sigma contamination
%             Sigma_outlier = k*this.generateA09cormatrix(p,0.95); % A09
%             replacement = tLocation;
             Sigma_outlier = tCorrelation;
            
            contamination = this.ContaminationType.generateContamination(contaminationDegree, p, replacement, tLocation, Sigma_outlier);		
            
            cindices = randperm(n, contaminationDegree);
            samples(cindices, :) = contamination;
        end
    end
    
end

