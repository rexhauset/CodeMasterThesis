function pathloss = largeScalePathloss(obj,rx,tx,isLOS)
%PATHLOSS 3GPP Model 
%   Calculates Pathloss according to formula given by 3GPP
%   Simplified effective height to be equal to 1
%   Might not be 100% accurate
    d_2d = calc2DDistance(rx,tx);
    d_3d = calcDistance(rx,tx);
    f_c = obj.getFrequency / 1e9;

    h_ut = rx.getHeight;
    h_bs = tx.getHeight;
    %Calculate Breakpoint            
    h_e = 1;
    d_BP = 4 * (h_bs-h_e) * (h_ut-h_e)*f_c*1e9/3e8;
    
    if strcmp(obj.situation,"UMa")
        pl_1 = 28+22*log10(d_3d)+20*log10(f_c);
        pl_2 = 28+40*log10(d_3d)+20*log10(f_c)-9*log10((d_BP^2+(h_ut-h_bs)^2)); 
        if d_2d<d_BP
            pl_LOS = pl_1;
        else
            pl_LOS = pl_2;
        end
        if isLOS
            pathloss = pl_LOS;
            sigma = 4;
        else
            pl_NLOS = 13.54+39.08*log10(d_3d)+20*log10(f_c)-0.6*(h_ut-1.5);
            pathloss = max(pl_LOS,pl_NLOS);
            sigma = 6;
        end
        
    elseif strcmp(obj.situation,"UMi")
        pl_1 = 32.4+21*log10(d_3d)+20*log10(f_c);
        pl_2 = 32.4+40*log10(d_3d)+20*log10(f_c)-9.5*log10((d_BP^2+(h_ut-h_bs)^2)); 
        if d_2d<d_BP
            pl_LOS = pl_1;
        else
            pl_LOS = pl_2;
        end
        if isLOS
            pathloss = pl_LOS;
            sigma = 4;
        else
            pl_NLOS = 22.4+35.3*log10(d_3d)+21.3*log10(f_c)-0.3*(rx.getHeight-1.5);
            pathloss = max(pl_LOS,pl_NLOS);
            sigma = 0.782;
        end
    else
        throw(MException("MYEX:badSituation","Can't calculate PL for situation "+situation))
    end
%     pathloss = pathloss + normrnd(0,sigma); NOTE: pathloss handled in
%     smallScale Fading
end