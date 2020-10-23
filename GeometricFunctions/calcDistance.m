function distance = calcDistance(rx,tx)
%CALCDISTANEC Summary of this function goes here
%   Detailed explanation goes here
    distance = sqrt(sum((rx.getPosition-tx.getPosition).^2));
end

