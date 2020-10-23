classdef Receiver < TxRx
    %Class to represent a transmitter for all models
    
    properties
    end
    
    methods
        function obj = Receiver(x,y,z)
            obj@TxRx(x,y,z)
        end
    end
end
