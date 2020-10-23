function params = smallScaleParams(obj, rx,tx,isLOS)
%SMALLSCALEPARAMS Calculates the random parameters refering to clusters
%   Output is a map with all the parameters needed later on, as given by table 7.5-6 
%   Map contains:
%{
    Delayspread (DS), Azimuth spread of departure (ASD), 
    Azimuth spread of Arrival (ASA), zenith spread of departure (ZSD),
    Zenith spread of arrival (ZSA), Shadow Fading (SF) and
%}
    %Things To do: copy entire table for each of the 4 cases
    %Draw DS, ASD, ASA, ZSD, ZSA, SF, K (K only if needed) and put into map
    % All the other params into map
    fc = obj.getFrequency/1e9;
    d2d = calc2DDistance(rx,tx);
    h_ut = rx.getHeight();
    h_bs = tx.getHeight();
    
    if strcmp(obj.situation,"UMa")
        if isLOS
            muDS = -6.955-0.0963*log10(fc);
            sigDS = 0.66;
            
            muASD = 1.06+0.1114*log10(fc);
            sigASD = 0.28;
            
            muASA = 1.81;
            sigASA = 0.2;
            
            muZSA = 0.95;
            sigZSA = 0.16;
            
            muZSD = max(-0.5, -2.1*(d2d/1000)-0.01*(h_ut-1.5)+0.75);
            sigZSD = 0.4;
            offsetZSD = 0;
            
            muSF = 0;
            sigSF = 4;
            
            muK = 9;
            sigK = 3.5;
            
            ccASDvDS = 0.4;
            ccASAvDS = 0.8;
            ccASAvSF = -0.5;
            ccASDvSF = -0.5;
            ccDSvSF = -0.4;
            ccASDvASA = 0;
            ccASDvK = 0;
            ccASAvK = -0.2;
            ccDSvK = -0.4;
            ccSFvK = 0;
            
            ccZSDvSF = 0;
            ccZSAvSF = -0.8;
            ccZSDvK = 0;
            ccZSAvK = 0;
            ccZSDvDS = -0.2;
            ccZSAvDS = 0;
            ccZSDvASD = 0.5;
            ccZSAvASD= 0;
            ccZSDvASA = -0.3;
            ccZSAvASA = 0.4;
            ccZSDvZSA = 0;
            
            delayScaling = 2.5;
            
            muXPR = 8;
            sigXPR = 4;
            
            N_clusters = 12;
            M_rpcluster = 20;
            
            clusterDS = max(0.25,6.5622-3.4084*log10(fc));
            clusterASD = 5;
            clusterASA = 11;
            clusterZSA = 7;
            
            clusterShadowingStd = 3;
            
            corrDist = containers.Map;
            corrDist("DS") = 30;
            corrDist("ASD") = 18;
            corrDist("ASA") = 15;
            corrDist("SF") = 37;
            corrDist("K") = 12;
            corrDist("ZSA") = 15;
            corrDist("ZSD") = 15;          
        else %UMa NLOS
            muDS = -6.28-0.204*log10(fc);
            sigDS = 0.39;
            
            muASD = 1.5-0.1114*log10(fc);
            sigASD = 0.28;
            
            muASA = 2.08-0.27*log10(fc);
            sigASA = 0.11;
            
            muZSA = -0.3236*log10(fc)+1.512;
            sigZSA = 0.16;
            
            muZSD = max(-0.5, -2.1*(d2d/1000)-0.01*(h_ut-1.5)+0.9);
            sigZSD = 0.49;
            offsetZSD = 7.66*log10(fc)-5.96- ...
                10^((0.208*log10(fc)-0.782)*log10(max(25,d2d))- ...
                0.13*log10(fc)+2.03-0.07*(h_ut-1.5));
            
            muSF = 0;
            sigSF = 6;
            
            muK = NaN;
            sigK = NaN;
            
            ccASDvDS = 0.4;
            ccASAvDS = 0.6;
            ccASAvSF = 0;
            ccASDvSF = -0.6;
            ccDSvSF = -0.4;
            ccASDvASA = 0.4;
            ccASDvK = NaN;
            ccASAvK = NaN;
            ccDSvK = NaN;
            ccSFvK = NaN;
            
            ccZSDvSF = 0;
            ccZSAvSF = -0.4;
            ccZSDvK = NaN;
            ccZSAvK = NaN;
            ccZSDvDS = -0.5;
            ccZSAvDS = 0;
            ccZSDvASD = 0.5;
            ccZSAvASD= -0.1;
            ccZSDvASA = 0;
            ccZSAvASA = 0;
            ccZSDvZSA = 0;
            
            delayScaling = 2.3;
            
            muXPR = 7;
            sigXPR = 3;
            
            N_clusters = 20;
            M_rpcluster = 20;
            
            clusterDS = max(0.25,6.5622-3.4084*log10(fc));
            clusterASD = 2;
            clusterASA = 15;
            clusterZSA = 7;
            
            clusterShadowingStd = 3;
            
            corrDist = containers.Map;
            corrDist("DS") = 40;
            corrDist("ASD") = 50;
            corrDist("ASA") = 50;
            corrDist("SF") = 50;
            corrDist("K") = NaN;
            corrDist("ZSA") = 50;
            corrDist("ZSD") = 50;
        end
        
    elseif strcmp(obj.situation,"UMi")
        if isLOS
            muDS = -0.24*log10(1+fc)-7.14;
            sigDS = 0.38;
            
            muASD = -0.05*log10(1+fc)+1.21;
            sigASD = 0.41;
            
            muASA = -0.08*log10(1+fc)+1.73;
            sigASA = 0.014*log10(1+fc)+0.28;
            
            muZSA = -0.1*log10(1+fc)+0.73;
            sigZSA = -0.04*log10(1+fc)+0.34;
            
            muZSD = max(-0.21, -14.8*(d2d/1000)-0.01*abs(h_ut-h_bs)+0.83);
            sigZSD = 0.35;
            offsetZSD = 0;
            
            muSF = 0;
            sigSF = 4;
            
            muK = 9;
            sigK = 5;
            
            ccASDvDS = 0.5;
            ccASAvDS = 0.8;
            ccASAvSF = -0.4;
            ccASDvSF = -0.5;
            ccDSvSF = -0.4;
            ccASDvASA = 0.4;
            ccASDvK = -0.2;
            ccASAvK = -0.3;
            ccDSvK = -0.7;
            ccSFvK = 0.5;
            
            ccZSDvSF = 0;
            ccZSAvSF = 0;
            ccZSDvK = 0;
            ccZSAvK = 0;
            ccZSDvDS = 0;
            ccZSAvDS = 0.2;
            ccZSDvASD = 0.5;
            ccZSAvASD= 0.3;
            ccZSDvASA = 0;
            ccZSAvASA = 0;
            ccZSDvZSA = 0;
            
            delayScaling = 3;
            
            muXPR = 9;
            sigXPR = 3;
            
            N_clusters = 12;
            M_rpcluster = 20;
            
            clusterDS = 5;
            clusterASD = 3;
            clusterASA = 17;
            clusterZSA = 7;
            
            clusterShadowingStd = 3;
            
            corrDist = containers.Map;
            corrDist("DS") = 7;
            corrDist("ASD") = 8;
            corrDist("ASA") = 8;
            corrDist("SF") = 10;
            corrDist("K") = 15;
            corrDist("ZSA") = 12;
            corrDist("ZSD") = 12;          
        else %UMi NLOS
            muDS = -0.24*log10(1+fc)-6.83;
            sigDS = 0.16*log10(1+fc)+0.28;
            
            muASD = -0.23*log10(1+fc)+1.53;
            sigASD = 0.11*log10(1+fc)+0.33;
            
            muASA = -0.08*log10(1+fc)+1.81;
            sigASA = 0.05*log10(1+fc)+0.3;
            
            muZSA = -0.04*log10(1+fc)+0.92;
            sigZSA = -0.07*log10(1+fc)+0.41;
            
            muZSD = max(-0.5, -3.1*(d2d/1000)-0.01*max(h_ut-h_bs,0)+0.2);
            sigZSD = 0.35;
            offsetZSD = -10^(-1.5*log10(max(10,d2d))+3.3);
            
            muSF = 0;
            sigSF = 7.82;
            
            muK = NaN;
            sigK = NaN;
            
            ccASDvDS = 0;
            ccASAvDS = 0.4;
            ccASAvSF = -0.4;
            ccASDvSF = 0;
            ccDSvSF = -0.7;
            ccASDvASA = 0;
            ccASDvK = NaN;
            ccASAvK = NaN;
            ccDSvK = NaN;
            ccSFvK = NaN;
            
            ccZSDvSF = 0;
            ccZSAvSF = 0;
            ccZSDvK = NaN;
            ccZSAvK = NaN;
            ccZSDvDS = -0.5;
            ccZSAvDS = 0;
            ccZSDvASD = 0.5;
            ccZSAvASD= 0.5;
            ccZSDvASA = 0;
            ccZSAvASA = 0.2;
            ccZSDvZSA = 0;
            
            delayScaling = 2.1;
            
            muXPR = 8;
            sigXPR = 3;
            
            N_clusters = 19;
            M_rpcluster = 20;
            
            clusterDS = 11;
            clusterASD = 10;
            clusterASA = 22;
            clusterZSA = 7;
            
            clusterShadowingStd = 3;
            
            corrDist = containers.Map;
            corrDist("DS") = 10;
            corrDist("ASD") = 10;
            corrDist("ASA") = 9;
            corrDist("SF") = 13;
            corrDist("K") = NaN;
            corrDist("ZSA") = 10;
            corrDist("ZSD") = 10;
        end
    else
        throw(MException("MYEX:badSituation","Can't calculate cluster parameters for situation "+situation))
    end
    
    % Draw random variables 
    mus = [muSF muK muDS muASD muASA muZSD muZSA];
    sigs = [sigSF sigK sigDS sigASD sigASA sigZSD sigZSA];
%     sig = diag([sigSF sigK 10^sigDS 10^sigASD 10^sigASA 10^sigZSD 10^sigZSA].^2);
    ccMatrix = diag([sigSF sigK sigDS sigASD sigASA sigZSD sigZSA].^2);
    ccMatrix = diag([1 1 1 1 1 1 1]);
       
    ccMatrix(1,2:end) = [ccSFvK ccDSvSF ccASDvSF ccASAvSF ccZSDvSF ccZSAvSF];
    ccMatrix(2,3:end) = [ccDSvK ccASDvK ccASAvK ccZSDvK ccZSAvK];
    ccMatrix(3,4:end) = [ccASDvDS ccASAvDS ccZSDvDS ccZSAvDS];
    ccMatrix(4,5:end) = [ccASDvASA ccZSDvASD ccZSAvASD];
    ccMatrix(5,6:end) = [ccZSDvASA ccZSAvASA];
    ccMatrix(6,7:end) = ccZSDvZSA;
    
    ccMatrix = (ccMatrix+ccMatrix') - eye(size(ccMatrix,1)).*diag(ccMatrix);    
    
% Second approach, closer to Winner model, but missing the Grid and filter
% part..: Cholesky decomposition of covariance matrix
    
    covMatrix = zeros(7,7);
    for i = 1:7
        for j = 1:7
            covMatrix(i,j) = ccMatrix(i,j) * sigs(i) * sigs(j);
        end
    end
    
    if isLOS
        vals = mvnrnd(zeros(7,1),[1 1 1 1 1 1 1]);
        vals = vals * chol(covMatrix);
        vals = vals + mus;
        k_offset = 1;
        K = vals(2);
    else
        % K not applicable
        covMatrix(2,:)=[];
        covMatrix(:,2)=[];
        mus(2)=[];
        sigs(2)=[];
        vals = mvnrnd(zeros(6,1),[1 1 1 1 1 1]);
        vals = vals * chol(covMatrix);
        vals = vals + mus;

        K = NaN;
        k_offset = 0;
    end
    
    
%    First approach, probably wront
%     if isLOS
%         vals = mvnrnd(mus, ccMatrix);
%         K = vals(2);
%         k_offset = 1;
%     else
%         % K not applicable
%         ccMatrix(2,:)=[];
%         ccMatrix(:,2)=[];
%         mus(2)=[];
%         vals = mvnrnd(mus, ccMatrix);
%         K = NaN;
%         k_offset = 0;
%     end
    
%     SF = vals(1);
%     DS = vals(2+k_offset);
%     ASD = vals(3+k_offset);
%     ASA = vals(4+k_offset);
%     ZSD = vals(5+k_offset);
%     ZSA = vals(6+k_offset);
    
    SF = vals(1);
    DS = 10^vals(2+k_offset);
    ASD = 10^vals(3+k_offset);
    ASA = 10^vals(4+k_offset);
    ZSD = 10^vals(5+k_offset);
    ZSA = 10^vals(6+k_offset);
    
    
    % Limit some of the values as specified by the standard
    ASA = min(ASA, 104);
    ASD = min(ASD, 104);
    ZSA = min(ZSA,52);
    ZSD = min(ZSD, 52);
    
    % Put all values into a map object
    
    params = containers.Map;
    params('DS') = DS;
    params('ASD') = ASD;
    params('ASA') = ASA;
    params('ZSA') = ZSA;
    params('ZSD') = ZSD;
    params('SF') = SF;
    params('K') = K;
 
    params('muZSD') = muZSD;
    params('offsetZSD') = offsetZSD;
    
    params('delayScaling') = delayScaling;
    params('muXPR') = muXPR;
    params('sigXPR') = sigXPR;
    params('N_clusters') = N_clusters;
    params('M_rpcluster') = M_rpcluster;
    
    params('clusterDS') = clusterDS;
    params('clusterASD') = clusterASD;
    params('clusterASA') = clusterASA;
    params('clusterZSA') = clusterZSA;
    
    params('clusterShadowingStd') = clusterShadowingStd;

    params('corrDist') = corrDist;    
end

