function labelPowerSpectrum(steps,nPower,minValue,dt,N)
% labelPowerSpectrum plots and labels the normalized power spectrum returned by plotPowerSpectrum
%
% SYNOPSIS    labelPowerSpectrum(steps,nPower,minValue,dt)
%
% INPUT       signal   : vector of normalized positions as returned by plotPowerSpectrum
%             s        : power spectrum vector as returned by plotPowerSpectrum
%             minValue : value [0..1] below which peak are not labeled
%             dt       : sampling interval (in seconds)
%
% DEPENDENCES labelPowerSpectrum uses { }
%             labelPowerSpectrum is used by {}
%
% Aaron Ponti, 2002

% Checki input parameters
if nargin==3
    time=0;
elseif nargin==4
    N=2*length(nPower);
    time=1;
elseif nargin==5
    time=1;
else
    error('Wrong number of input parameters');
end
    
% Plot the power spectrum
figure;
plot(steps,nPower,'k-');
% stem(steps,nPower);
hold on;

% Read first value
pos=find(nPower==max(nPower));
value=nPower(pos);

% Go through peaks and label them if their power is > that minValue
while value>minValue
    
    if time==1
        if pos~=1 % We don't label the DC
            period=N/pos*dt;
            text(steps(pos),value,[num2str(round(period)),'s']);
        end
    else
        if pos~=1 % We don't label the DC
            text(steps(pos),value,num2str(pos));
        end
    end
    nPower(pos)=0;

    % Read next value
    pos=find(nPower==max(nPower));
    value=nPower(pos);

end

hold off;