function coefficients = H_usnNLOS(clusterPowers, F_rx, F_tx, rx, tx, kappa, phases, rx, tx, fc, phi_nmAOA, theta_nmZOA, phi_nmAOD, theta_nmZOD)
%H_USNNLOS Calculates channel coefficients for clusters 3 ... N
%   Formula 7.5-22
    
    M_numberRays = size(kappa,2);
    lambda = 3e8/fc;
    
    rrx = [sind(theta_nmZOA) * cosd(phi_nmAOA);
        sind(theta_nmZOA) * sind(phi_nmAOA);
        cosd(theta_nmZOA)];
    drx = rx.getPosition();
    
    rtx = [sind(theta_nmZOD) * cosd(phi_nmAOD);
        sind(theta_nmZOD) * sind(phi_nmAOD);
        cosd(theta_nmZOD)];
    dtx = tx.getPosition();
    
    polarizationMatrix = exp(j*phases);
    for i = 1:size(clusterPowers,1);
        for j = 1:M_numberRays
            polarizationMatrix(i,j) = polarizationMatrix(i,j) .* [1 sqrt(kappa(i,j)^-1;sqrt(kappa(i,j)^-1 1];
    
    Husn =  sqrt(clusterPowers/M_numberRays) * F_rx * polarizationMatrix * F_tx * ...
        exp(j*2*pi*rrx*drx/lambda) * exp(j*2*pi*rtx*dtx/lambda);
end

