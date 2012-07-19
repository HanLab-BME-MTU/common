function [PersTime,newBlock] = findingProtRetrTime(block,TS,deltaT)
%This is the core function to calculate the protrusion/retrection time
%It is part of a set of function (see below)
%
% Usage: [PersTime,newBlock] = findingProtRetrTime(block,TS,deltaT)
%
% Input: block  - cell array with the time blocks - (see getPersistenceTime )
%        TS     - Time series
%        deltaT - Time between points (Sampling rate)
%
%Output: PersTime - vector with the persistence time for each block
%        newBlock - Two blocks are fused when they are separeted by one point and this point does not cross the x-axis.
%
% See also: getPersistenceTime, getEdgeMotionPersistence
%
%Marco Vilela, 2012

%% Parsing the input ******************************************************
ip = inputParser;
ip.addRequired('block',@iscell);
ip.addRequired('TS',@isvector);
ip.addRequired('deltaT',@isscalar);
%**************************************************************************
    nB     = length(block);
    nPoint = length(TS);

    %Fusing blocks that are 1 point apart - The point has to be tested
    newBlock{1} = block{1};
    cc = 1;
    for iB = 2:nB
        if ( block{iB}(1) - newBlock{cc}(end) == 2 ) & TS( block{iB}(1) - 1 ) > 0%Testing gap
            newBlock{cc} = [newBlock{cc};block{iB}(1)-1;block{iB}];
        else
            cc = cc + 1;
            newBlock{cc} = block{iB};
        end
    end

    PersTime  = NaN(cc,1);

    for iB = 1:cc

        if newBlock{iB}(1) == 1

            deltaTr = rightBorder(TS,newBlock{iB},deltaT,nPoint);
            deltaTl = 0;

        elseif newBlock{iB}(end) == nPoint

            deltaTl = leftBorder(TS,newBlock{iB},deltaT);
            deltaTr = 0;

        else

            deltaTr = rightBorder(TS,newBlock{iB},deltaT,nPoint);
            deltaTl = leftBorder(TS,newBlock{iB},deltaT);

        end

        PersTime(iB)  = (length(newBlock{iB})-1)*deltaT + deltaTl + deltaTr;


        clear deltaTl;clear deltaTr;
    end

end%End of main function

function [deltaTr] = rightBorder(TS,newBlock,deltaT,nPoint)

        endPoint = find( TS(newBlock(end)+1:end) < 0 );%finding the zero-crossing point
        
        if ~isempty(endPoint)
            rZeroC2  = endPoint(1) + newBlock(end);
            rZeroC1  = rZeroC2  - 1;
            deltaTr  = deltaT*TS(rZeroC1)/abs( diff( TS(rZeroC1:rZeroC2) ) );
        else
            rZeroC2  = nPoint;
            rZeroC1  = newBlock(end);
            deltaTr  = deltaT*(rZeroC2 - rZeroC1);
        end
       deltaTr = deltaTr + deltaT*(rZeroC1 - newBlock(end)); 
end

function [deltaTl] = leftBorder(TS,newBlock,deltaT)
    
    ftPoint  = find( TS(1:newBlock(1)) < 0 );%finding the zero-crossing point
    
    if ~isempty(ftPoint)
        lZeroC1  = ftPoint(end) ;
        lZeroC2  = lZeroC1  + 1;
        deltaTl  = deltaT*TS(lZeroC2)/abs( diff( TS(lZeroC1:lZeroC2) ) );
    else
        lZeroC1  = 1;
        lZeroC2  = newBlock(end);
        deltaTl  = deltaT*(lZeroC2 - lZeroC1);
    end
    deltaTl = deltaTl + deltaT*(newBlock(1) - lZeroC2);
end