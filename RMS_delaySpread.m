function rmsDS = RMS_delaySpread(powers, times)
%Function to calulate RMS Delay Spread
% Expects two 1D vectors. The first containing powers, the 2nd the
% corresponding delays.
    
    assert(size(powers,1)==size(times,1),"Sizes of powers and times must be identical)")
    assert(size(powers,2)==size(times,2),"Sizes of powers and times must be identical)")
    assert(size(powers,2)==1, "Arrays should be 1D")
    
    excessDelay_1 = (powers' * times) / sum(powers);
    excessDelay_2 = (powers' * times.^2) / sum(powers);
    rmsDS =  sqrt(excessDelay_2 - excessDelay_1^2);
end

