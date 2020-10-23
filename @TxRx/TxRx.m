classdef TxRx < handle
    %TXRX Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        position
    end
    
    methods
        function obj = TxRx(x,y,z)
            %TXRX Construct an instance of this class
            obj.position = [x,y,z];
        end
        
        function obj = movePosition(obj,dx,dy,dz)
            obj.position = obj.position + [dx,dy,dz];
        end
        
        function pos = getPosition(obj)
            pos = obj.position;
        end
        
        function setPosition(obj,x,y,z)
            obj.position = [x,y,z];
        end
        
        function h=getHeight(obj)
            h=obj.position(3);
        end
        
        radiatedPower = antennaDirectionalPattern(obj, phi, theta)

    end
end

