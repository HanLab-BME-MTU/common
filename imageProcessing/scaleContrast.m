function out = scaleContrast(in, varargin)
switch(nargin)
    case 1
        maxI = max(in(:));
        minI = min(in(:));
        v = 255;
    case 2
        dynamicRange = varargin{1};
        minI = dynamicRange(1);
        maxI = dynamicRange(2);
        v = 255;
    case 3
    otherwise
        disp('Unsupported number of arguments');
end;
if ((maxI-minI) ~= 0)
    out = (in-minI)/(maxI-minI)*v;
else
    out = zeros(size(in));
end;