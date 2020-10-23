function [AOD, ZOD, AOA, ZOA] = anglesBetweenRxTx(rx, tx)
%ANGLESBETWEENRXTX Calculates azimuths and Zeniths of arrival and departure
%   To be reviewed for all cases ...

%   Zeniths:
    d2d = calc2DDistance(rx,tx);
    dz = abs(tx.getHeight-rx.getHeight);
    
    ZOD_inner = rad2deg(atan(d2d/dz));
    ZOA = ZOD_inner;
    ZOD = 180 - ZOD_inner;
    
    rxPos = rx.getPosition;
    txPos = tx.getPosition;
    
%     Azimuths
    dx = rxPos(1)-txPos(1);
    dy = rxPos(2)-txPos(2);
    
    X = [1 0 0];
    v = rxPos - txPos;
    
    CosTheta = max(min(dot(X,-v)/(norm(X)*norm(v)),1),-1);
    AOA = real(acosd(CosTheta));
    AOD = 180-AOA;
    
    if dy < 0
        AOD = -AOD;
    elseif dy > 0
        AOA = -AOA;
    end
%     if dy == 0
%         if rxPos(1) < rxPos(2)
%             AOD = 180;
%         else
%             AOD = 0;
%         end
%     else
%         AOD = rad2deg(atan(dx/dy));
%         AOA = 180 - AOD;
%     end
       
end