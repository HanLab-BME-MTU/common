function nT = numTimePoints(obj)
% get the number of time points available
        M = vertcat(obj.seqOfEvents);
        nT = max(M(:,1));
end
