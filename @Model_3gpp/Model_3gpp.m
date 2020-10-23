classdef Model_3gpp < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        situation % can be UMi or UMa
        frequency
        bandwidth
        pathloss
    end
    
    methods
        function obj = Model_3gpp(situation)
            if nargin == 0
                obj.situation = "UMi";
            else
                obj.situation = situation;
            end
            obj.frequency = 28e9;
            obj.bandwidth = 100e6;
        end
        pathloss = largeScalePathloss(obj,rx,tx,isLOS)
        [isLOS,pLOS] = rollLOS(obj,rx,tx)
        params = smallScaleParams(obj, rx,tx,isLOS)
        fading = applyModel(obj,rx,tx,isLOS)
        function freq = getFrequency(obj)
            freq = obj.frequency;
        end
        
    end
end

