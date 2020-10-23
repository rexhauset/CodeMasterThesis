function results = applyModel(obj,rx,tx,isLOSOverride)
%SMALLSCALEFADING Entire small scale fading happens here
%   Implementation of Steps 1 to 11 in standard


% Step 1: Environment setup
    % Determine LOS angles: AOD ZOD AOA ZOA
    [AOD, ZOD, AOA, ZOA] = anglesBetweenRxTx(rx, tx);
    
    % Field Patters set in TxRx Class
    
    % Angles of antenna and Base station !? TODO
    
    % Speed and direction of UT
%     v_ut = 0;
    
% Step 2: LOS or NLOS?
    isLOS = obj.rollLOS(tx,rx); 
    if exist('isLOSOverride','var')==1
        isLOS = isLOSOverride;
    end
% Step 3: Pathloss
    obj.pathloss = obj.largeScalePathloss(tx,rx,isLOS);

    
% Step 4: Gemerate large scale parameters
    parameters = obj.smallScaleParams(rx,tx,isLOS);

% SMALL SCALE PARAMETERS
% Step 5: Cluster delays
    N_clusters = parameters('N_clusters');
    delayScalingParamter = parameters('delayScaling');
    delaySpread = parameters('DS');
    
    clusterDelays = - delayScalingParamter * delaySpread * log(rand(N_clusters,1));
    clusterDelays = clusterDelays - min(clusterDelays);
    clusterDelays = sort(clusterDelays);
    
    if isLOS
        K = parameters('K');
        C_tau = 0.7705 - 0.0433 * K + 0.0002 * K^2 + 0.000017 * K^3;
        clusterDelays_LOS = clusterDelays / C_tau;
    end
    
% Step 6: Cluster Powers 
    perClusterShadowingVar = parameters('clusterShadowingStd');
    perClusterShadowing = normrnd(0,perClusterShadowingVar^2);
    
    clusterPowers = exp(-clusterDelays * (delayScalingParamter-1) / (delayScalingParamter*delaySpread)) * 10^(-perClusterShadowing/10);
    
    clusterPowers_n = clusterPowers / sum(clusterPowers);
    clusterPowers75_22 = clusterPowers_n; %Copy of cluster powers to be used in Eq. 7.5-22
    
    if isLOS
        K_R = 10^(K / 10);
        P1LOS = K_R / (K_R + 1);
        
        clusterPowers = 1 / (K_R + 1) * clusterPowers / sum(clusterPowers);
        clusterPowers(1) = clusterPowers(1) + P1LOS;
    else
        clusterPowers = clusterPowers_n;
    end
    
    M_raysPercluster = parameters('M_rpcluster');
    rayPowers = clusterPowers / M_raysPercluster;
    rayPowers = diag(rayPowers) * ones(N_clusters, M_raysPercluster);
    
    % Delete cluster with less than -25dB compared to max
    maxPower = max(clusterPowers);
    toRetain = maxPower / clusterPowers <= 10^2.5;
    
    clusterPowers = clusterPowers(toRetain);
    clusterPowers75_22 = clusterPowers75_22(toRetain);
    rayPowers = rayPowers(toRetain,:);
        
    N_clusters = size(clusterPowers,1);
    
    % Set cluster delays to incorporate LOS changes
    if isLOS
        clusterDelays = clusterDelays_LOS;
    end
    
 %Step 7: Arrival angles and departure angles for clusters
    % TODO: C_phi for angles is weird, only given for some values. So far,
    % simply used linear approximation.
    
    % Starting with Azimuths
    % Declarations for arrival and departure:
    C_phis = [0 0 0 0.779 0.86 0.915 0.965 1.018 1.058 1.09 1.123 1.146 1.168 1.190 1.211 1.226 1.243 1.257 1.273 1.289];
    assert(N_clusters>3,"Not sufficient data available to model for less than 3 clusters.")
    C_phi = C_phis(N_clusters);
    
    rayOffsets = [0.0447 -0.0447 0.1413 -0.1413 0.2492 -0.2492 ...
        0.3715 -0.3715 0.5129 -0.5129 0.6797 -0.6797 ...
        0.8844 -0.8844 1.1481 -1.1481 1.5195 -1.5195 2.1551 -2.1551];
    ASA = parameters('ASA');
    ASD = parameters('ASD');
    clusterASA = parameters('clusterASA');
    clusterASD = parameters('clusterASD');
    
    
    if isLOS
        C_phi = C_phi * (1.1035 - 0.058 * K - 0.002 * K^2 + 0.0001 * K^3);
    end
    
    %%% AOA %%%
    phi_nAOA = 2 * (ASA / 1.4) * sqrt(-log(clusterPowers/max(clusterPowers))) / C_phi;
    
    
    Yn = normrnd(0,(ASA/7)^2,N_clusters,1);
    Xn = (randi(2,N_clusters,1)*2-3);
    
    if isLOS
        phi_nAOA = phi_nAOA .* Xn + Yn + AOA - Xn(1) .* phi_nAOA(1) - Yn(1);
    else
        phi_nAOA = phi_nAOA .* Xn + Yn + AOA;
    end
    phi_nAOA = phi_nAOA - fix(phi_nAOA./360)*360;
    phi_nAOA(abs(phi_nAOA) > 180) = - sign(phi_nAOA(abs(phi_nAOA) > 180)) .* (360 - abs(phi_nAOA(abs(phi_nAOA)>180)));
    
    phi_nmAOA = phi_nAOA + ones(N_clusters,1) * rayOffsets * clusterASA;
    phi_nmAOA = phi_nmAOA - fix(phi_nmAOA./360)*360;
    phi_nmAOA(abs(phi_nmAOA) > 180) = - sign(phi_nmAOA(abs(phi_nmAOA) > 180)) .* (360 - abs(phi_nmAOA(abs(phi_nmAOA) > 180)));
    %%% AOD %%%
    phi_nAOD = 2 * (ASD / 1.4) * sqrt(-log(clusterPowers/max(clusterPowers))) / C_phi;
    
    
    Yn = normrnd(0,(ASD/7)^2,N_clusters,1);
    Xn = (randi(2,N_clusters,1)*2-3);
    
    if isLOS
        phi_nAOD = phi_nAOD .* Xn + Yn + AOD - Xn(1) .* phi_nAOD(1) - Yn(1);
    else
        phi_nAOD = phi_nAOD .* Xn + Yn + AOD;
    end
    
    phi_nAOD = phi_nAOD - fix(phi_nAOD./360)*360;
    phi_nAOD(abs(phi_nAOD) > 180) = - sign(phi_nAOD(abs(phi_nAOD) > 180)) .* (360 - abs(phi_nAOD(abs(phi_nAOD) > 180)));

    
    phi_nmAOD = phi_nAOD + ones(N_clusters,1) * rayOffsets * clusterASD;
    phi_nmAOD = phi_nmAOD - fix(phi_nmAOD./180)*180;
    phi_nmAOD(abs(phi_nmAOD) > 180) = - sign(phi_nmAOD(abs(phi_nmAOD) > 180)) .* (360 - abs(phi_nmAOD(abs(phi_nmAOD) > 180)));

    % Zeniths
    ZSA = parameters('ZSA');
    ZSD = parameters('ZSD');
    clusterZSA = parameters('clusterZSA');
    offsetZSD = parameters('offsetZSD');
    muZSD = parameters('muZSD');
    
    % Values approximated, linear regression and comparison with Azimuth
    % factors
    C_phis = [0 0 0 0.579 0.66 0.715 0.765 0.889 0.92 0.957 1.031 1.104 1.105 1.107 1.1088 1.125 1.142 1.16 1.184 1.178];
%     assert(N_clusters>3,"Not sufficient data available to model for less than 3 clusters.")
    C_phi = C_phis(N_clusters);

    if isLOS
        C_phi = C_phi * (1.3086 + 0.0339 * K - 0.0077 * K^2 + 0.0002 * K^3);
    end
        
    %%% ZOA %%%
    theta_nZOA = -ZSA * log(clusterPowers/max(clusterPowers)) / C_phi;
    
    Yn = normrnd(0,(ZSA/7)^2,N_clusters,1);
    Xn = (randi(2,N_clusters,1)*2-3);
    
    if isLOS
        theta_nZOA = Xn .* theta_nZOA + Yn + ZOA - Xn(1) * theta_nZOA(1) - Yn(1);
    else
        theta_nZOA = Xn .* theta_nZOA + Yn + ZOA;
    end
    
    theta_nmZOA = theta_nZOA + ones(N_clusters,1) * rayOffsets * clusterZSA;
    theta_nmZOA(theta_nmZOA>180) = 360 - theta_nmZOA(theta_nmZOA>180);
    
    %%% ZOD %%%
    theta_nZOD = -ZSA * log(clusterPowers/max(clusterPowers)) / C_phi;
    
    Yn = normrnd(0,(ZSD/7)^2,N_clusters,1);
    Xn = (randi(2,N_clusters,1)*2-3);
    
    if isLOS
        theta_nZOD = Xn .* theta_nZOD + Yn + ZOD - Xn(1) * theta_nZOD(1) - Yn(1);
    else
        theta_nZOD = Xn .* theta_nZOD + Yn + ZOD + offsetZSD;
    end
    
    theta_nmZOD = theta_nZOD + ones(N_clusters,1) * rayOffsets * 3/8 * 10^muZSD;
    theta_nmZOD(theta_nmZOD>180) = 360 - theta_nmZOD(theta_nmZOD>180);
    
 % Step 8: Coupling of angles
 % Initial idea: shuffle all nm Matrices columns and assume: same column = coupled.
 % Has to be revisited with Step 11 maybe.
    cols = M_raysPercluster;
    indices = randperm(cols);
    theta_nmZOD = theta_nmZOD(:,indices);
    indices = randperm(cols);
    phi_nmAOA = phi_nmAOA(:,indices);
    indices = randperm(cols);
    theta_nmZOA = theta_nmZOA(:,indices);
    indices = randperm(cols);
    phi_nmAOD = phi_nmAOD(:,indices);
    
 % Step 9: Generate Cross Polarization Power Ratios (XPR, kappa_nm)
    muXPR = parameters('muXPR');
    sigXPR = parameters('sigXPR');
    
    kappa_nm = 10.^(normrnd(muXPR, sigXPR^2, N_clusters, M_raysPercluster)/10);
    
% Coefficient Generation
 % Step 10: Draw initial random phases
    % Convention:
    % 1st dimension: cluster, 2nd dimension: ray
    % 3rd: 1: theta theta, 2: theta phi, 3: phi theta 4: phi phi
    initialPhases_polar_nm = rand(N_clusters, M_raysPercluster, 4)*2*pi - pi;
    
 % Step 11: Generate Channel coefficients
    AOAs = AOA - phi_nmAOA;
    ZOAs = 90 - (ZOA - theta_nmZOA);
    
    AODs = AOD - phi_nmAOD;
    ZODs = 90 - (ZOD - theta_nmZOD);
    
    polarization_angle = pi/4;
    
    F_rx = 10.^(rx.antennaDirectionalPattern(AOAs, ZOAs)/10);
    F_rx_theta = sqrt(F_rx) * cos(polarization_angle);
    F_rx_phi = sqrt(F_rx) * sin(polarization_angle);
    
    F_tx = 10.^(tx.antennaDirectionalPattern(AODs, ZODs)/10);
    F_tx_theta = sqrt(F_tx) * cos(polarization_angle);
    F_tx_phi = sqrt(F_tx) * sin(polarization_angle);

    lambda = 3e8/obj.getFrequency();
    
    rrx = cat(3, sind(theta_nmZOA) .* cosd(phi_nmAOA),...
        sind(theta_nmZOA) .* sind(phi_nmAOA),...
        cosd(theta_nmZOA));
    drx = rx.getPosition();
    
    rtx = cat(3, sind(theta_nmZOD) .* cosd(phi_nmAOD),...
        sind(theta_nmZOD) .* sind(phi_nmAOD),...
        cosd(theta_nmZOD));
    dtx = tx.getPosition();
    
    polarizationMatrix = exp(1i .* initialPhases_polar_nm);
    polarMatrix = zeros(N_clusters,M_raysPercluster,2,2);
    
    Husnm = zeros(N_clusters,M_raysPercluster);
    Husn = zeros(N_clusters);
    
    % Calculation for Clusters 3 ... N.
    for n = 1:N_clusters
        for m = 1:M_raysPercluster
            polarMatrix(n,m,:,:) = reshape(polarizationMatrix(n,m,:),[2,2]) .* [1, sqrt(kappa_nm(n,m)^-1);sqrt(kappa_nm(n,m)^-1), 1];
            Husnm(n,m) = sqrt(clusterPowers75_22(n)/M_raysPercluster) * [F_rx_theta(n,m) F_rx_phi(n,m)] * reshape(polarMatrix(n,m,:,:),2,2) * [F_tx_theta(n,m); F_tx_phi(n,m)] * ...
                exp(1i*2*pi*reshape(rrx(n,m,:),1,3)*drx'/lambda) * exp(1i*2*pi*reshape(rtx(n,m,:),1,3)*dtx'/lambda);
            
        end
        if n > 2
            Husn(n) = sum(Husnm(n,:));
        end
    end
    
    % Calculations for Clusters 1,2:
    % Three subclusters 1, 2, 3. Ray indices:
    R{1} = [1:8,19,20];
    R{2} = [9:12,17];
    R{3} = 13:16;
    
    % Delays of subclusters:
    clusterDS = parameters('clusterDS');
    taus = [0, 1.28, 2.56] * clusterDS;
    
    Husni = zeros(2,3);
    Hust = 0;
%     syms t
    for n = 1:2
        for i = 1:3
            Husni(n,i) = sum(Husnm(n,R{i})); 
%             Hust = Hust + size(R{i},2) / M_raysPercluster * Husni(n,i) * dirac(t-(clusterDelays(n)+taus(i)));
        end
    end

    for n = 3:N_clusters
%         Hust = Hust + Husn(n) * dirac(t-clusterDelays(n));
    end
    
    if isLOS
        F_rx = rx.antennaDirectionalPattern(90,0);
        F_tx = tx.antennaDirectionalPattern(90,0);
        
        F_rx_theta = sqrt(F_rx) * cos(polarization_angle);
        F_rx_phi = sqrt(F_rx) * sin(polarization_angle);
        
        F_tx_theta = sqrt(F_tx) * cos(polarization_angle);
        F_tx_phi = sqrt(F_tx) * sin(polarization_angle);
        
        rrx = [sind(ZOA) * cosd(AOA),...
        sind(ZOA) * sind(AOA),...
        cosd(ZOA)];
    
        rtx = [sind(ZOD) * cosd(AOD),...
        sind(ZOD) * sind(AOD),...
        cosd(ZOD)];
        
%         Hus1 = [F_rx_theta F_rx_phi] * [1 0;0 -1] * [F_tx_theta; F_tx_phi] * ...
%                 exp(-1i*2*pi*calcDistance(rx,tx)/lambda) * exp(1i*2*pi*rrx*drx'/lambda) * exp(1i*2*pi*rtx*dtx'/lambda);
        
%         Hust = sqrt(1/(K_R+1)) * Hust + sqrt(K_R/(K_R+1)) * Hus1 * dirac(t-clusterDelays(1));
        
    end
    
    results = containers.Map;
    angles = containers.Map;
    angles('AOA') = phi_nAOA;
    angles('AOD') = phi_nAOA;
    angles('ZOA') = phi_nAOA;
    angles('ZOD') = phi_nAOA;
    angles('phi_nAOA') = phi_nAOA;
    angles('phi_nAOD') = phi_nAOD;
    angles('phi_nmAOA') = phi_nmAOA;
    angles('phi_nmAOD') = phi_nmAOD;
    angles('theta_nZOA') = theta_nZOA;
    angles('theta_nZOD') = theta_nZOD;
    angles('theta_nmZOA') = theta_nmZOA;
    angles('theta_nmZOD') = theta_nmZOD;

    results('N_clusters') = N_clusters;
    results('clusterDelays') = clusterDelays;
    results('clusterPowers') = clusterPowers;
    results('channel') = Hust;
    results('angles') = angles;
    results('isLOS') = isLOS;
    results('parameters') = parameters;

    % Formula 7.5-22 
%     Husn =  sqrt(clusterPowers/M_raysPercluster) * F_rx * polarMatrix * F_tx * ...
