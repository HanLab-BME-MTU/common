function particleBehaviorVsEdgeMotion(protPerWindow,sptPropInWindow,prop2analyze,band2analyze) %#ok<INUSL>
%PARTICLEBEHAVIORVSEDGEMOTION looks for correlation between particle behavior and cell edge protrusion activity
%
%SYNOPSIS 
%
%INPUT  protPerWindow  : Average protrusion vector per window (from
%                        Hunter).
%       sptPropInWindow: Output of particleBehaviorInWindows.
%       prop2analyze   : Cell array indicating properties to analyze. The
%                        names must exactly match the field names in
%                        sptPropInWindow.
%                        Optional. Default: All properties.
%       band2analyze   : Vector indicating which bands to analyze. Bands
%                        arise from the division of a cell into windows and
%                        run parallel to the cell edge. Band #1 is at the
%                        cell edge.
%                        Optional. Default: 1.
%
%OUTPUT 
%
%
%Khuloud Jaqaman, September 2010

%% Input

if nargin < 3 || isempty(prop2analyze)
    prop2analyze = {'spDensity','f2fDisp','fracUnclass','fracConf',...
        'fracBrown','fracDir','diffCoef','confRad'};
end

if nargin < 4 || isempty(band2analyze)
    band2analyze = 1;
end

%calculate number of properties and number of bands
numProp2analyze = length(prop2analyze);
numBands = length(band2analyze);

%% Plot activity maps of individual properties

%protrusion vectors
figure
protNormVecMag = protPerWindow.averageNormalComponent;
imagesc(protNormVecMag);
caxis([-1 1]*max(abs(protNormVecMag(:))))
colorbar
title('Protrusion Vectors')
xlabel('Time (frames)')
ylabel('Position along edge (window #)')

%go over requested properties
cmap = [[0 0 0]; colormap];
for iProp = 1 : numProp2analyze
    
    %get property name
    propName = prop2analyze{iProp};
    
    %get current property
    eval(['propCurrent = sptPropInWindow.' propName '.values;']);
    
    switch propName
        
        case {'spDensity','f2fDisp','diffCoef','confRad'}
            
            %go over all bands
            for iBand = 1 : numBands
                propCurrentBand = squeeze(propCurrent(band2analyze(iBand),:,:)); %#ok<NODEF>
                valuesNoNaNs = propCurrentBand(:);
                valuesNoNaNs = valuesNoNaNs(~isnan(valuesNoNaNs));
                minValue = prctile(valuesNoNaNs,1);
                maxValue = prctile(valuesNoNaNs,99);
                %                 minValue = min(propCurrentBand(:));
                %                 maxValue = max(propCurrentBand(:));
                nanValue = minValue - (maxValue - minValue)/64;
                propCurrentBand(isnan(propCurrentBand)) = nanValue;
                figure
                imagesc(propCurrentBand);
                caxis([nanValue maxValue])
                colormap(cmap);
                colorbar
                title(propName)
                xlabel('Time (frames)')
                ylabel('Position along edge (window #)')
            end
            
        case {'fracUnclass','fracConf','fracBrown','fracDir'}
            
            %go over all bands
            for iBand = 1 : numBands
                propCurrentBand = squeeze(propCurrent(band2analyze(iBand),:,:)); %#ok<NODEF>
                propCurrentBand(isnan(propCurrentBand)) = -1/64;
                figure
                imagesc(propCurrentBand);
                caxis([-1/64 1])
                colormap(cmap);
                colorbar
                title(propName)
                xlabel('Time (frames)')
                ylabel('Position along edge (window #)')
            end
            
    end
    
end






