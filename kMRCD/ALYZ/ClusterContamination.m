classdef ClusterContamination < handle
    %UNTITLED4 Summary of this class goes here
    %   Detailed explanation goes here    
    
    methods
        function xx = generateContamination(~, m, p, r, tLoc, tCor)
			xx = mvnrnd((r + tLoc)', tCor, m);
        end
    end    
end
