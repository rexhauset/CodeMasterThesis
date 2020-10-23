classdef Transmitter < TxRx
    %Class to represent a transmitter for all models
    
    properties
    end
    
    methods
        function obj = Transmitter(x,y,z)
            obj@TxRx(x,y,z)
        end
    end
end

