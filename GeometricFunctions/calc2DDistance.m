function distance = calc2DDistance(rx,tx)
%CALCDISTANEC Summary of this function goes here
%   Detailed explanation goes here
    d = rx.getPosition-tx.getPosition;
    d = d(1:2);
    distance = sqrt(sum(d.^2));
end

