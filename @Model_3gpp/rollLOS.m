function [isLOS, pLOS] = rollLOS(obj,rx,tx)
% Function to calculate whether or not there is a Line-of-Sight
% According to the 3GPP model
    d = calc2DDistance(rx,tx);

    if d < 18
        isLOS = true;
        pLOS = 1;
        return
    end
      
    if strcmp(obj.situation,"UMi")
        pLOS = 18/d + exp(-d/36) * (1-18/d);
    elseif strcmp(obj.situation,"UMa")
        f1 = 18/d + exp(-d/63) * (1-18/d);
        rxHeight = rx.getHeight();
        if rxHeight <= 13
            pLOS = f1;
        else
            C = ((rxHeight-13)/10)^1.5;
            f2 = 1+C*5/4*exp(-d/150)*(d/100)^3;
            pLOS = f1 * f2;
        end
    else
        throw(MException("MYEX:badSituation","Can't calculate LOS for situation "+situation))
    end

    isLOS = rand < pLOS;        
end