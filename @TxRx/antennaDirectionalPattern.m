function radiatedPower = antennaDirectionalPattern(obj, phi, theta)
%ANTENNA RADIATION PATTERN
%   According to 3GPP standard, page 22
%   theta, phi angles at antenna.
%   theta e 0,180 (deg), phi e -180,180 (deg)
    
    N_angles = size(phi,1) * size(phi,2);

    assert(sum(sum(theta <= 180 & theta >= 0))/N_angles, "RADIATION PATTERN ERROR: Theta must be between 0 and 180 degrees.")
    assert(sum(sum(abs(phi) <= 180))/N_angles, "RADIATION PATTERN ERROR: Phi must be between -180 and 180 degrees.")

    theta_3db = 65;
    phi_3db = 65;
    A_max = 8;
    vert = -min(12 * (theta-90).^2 / theta_3db^2, 30);
    hori = -min(12 * (phi).^2 / phi_3db^2, 30);
    
    radiatedPower = -min(-(vert+hori),A_max);
end

