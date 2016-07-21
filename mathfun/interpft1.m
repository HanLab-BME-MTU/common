function [vq] = interpft1(varargin)
% INTERPFT1 Evaluate 1D Fourier series at off-grid points, xq
%
% INPUT
% x - (optional) interval that the periodic function is defined on
%     default: 1 size(v,1)+1
% v - known values at regular intervals corresponding to
%     linspace(x(1),x(2)-x(1),size(v,1))
% xq - abisssa to determine value of Fourier series at
% method - (optional) interpolation method, see interp1
%          default: pchip
% fineGridFactor - (optional) Number of grid points to expand to using interpft
%                  default: Depends on method
%                  6 for cubic, pchip, v5cubic
%                  3 for spline
%                  10 for everything else
%
% OUTPUT
% vq - values of Fourier series at xq

% Mark Kittisopikul
% July 20th, 2016

% TODO:
% 1. Implement Horner's method, see trigtech.horner
% 2. Implement Matrix Multiplication Transform (MMT), see Boyd

    [x,v,xq,method,fineGridFactor] = parseinputs(varargin{:});
    vft3 = interpft(v,size(v,1)*fineGridFactor);
    vft3 = [vft3(end-2:end,:); vft3; vft3(1:4,:)];
    if(~isempty(x))
        % Map indices from [x(1) x(2)) to [0 1)
        period = x(2) - x(1);
        xq = xq - x(1);
        xq = mod(xq,period);
        xq = xq./period;
        % Map indices from [0 1) to [4 size(v,1)*fineGridFactor+4)
        xq = xq.*(size(v,1)*fineGridFactor);
        xq = xq+4;
    end
    vq = interp1(vft3,xq,method);
end

function [x,v,xq,method,fineGridFactor] = parseinputs(varargin)
    x = [];
    method = [];
    fineGridFactor = [];
    switch(nargin)
        case 2
            % v xq

            [v, xq]  = varargin{:};
        case 3
            if(ischar(varargin{3}))
                % v xq method
                [v, xq, method] = varargin{:};
            else
                % x v xq
                [x, v, xq] = varargin{:};
            end
        case 4
            if(ischar(varargin{4}))
                % x v xq method
                [x,v,xq,method] = varargin{:};
            else
                % v xq method fineGridFactor
                [v,xq,method,fineGridFactor] = varargin{:};
            end
        case 5
            [x, v,xq,method,fineGridFactor] = varargin{:};
        otherwise
            error('interpft1:nargin','Incorrect number of inputs');
    end
    if(isrow(v))
        v = v';
    end
    if(isempty(x))
        x = [1 size(v,1)+1];
    else
        assert(numel(x) == 2,'interpft1:IncorrectX', ...
            'x must either be empty or have 2 elements.');
    end
    if(isempty(method))
        method = 'pchip';
    end
    if(isempty(fineGridFactor))
        switch(method)
%             case 'linear'
%                 fineGridFactor = 10;
%             case 'nearest'
%                 fineGridFactor = 10;
%             case 'next'
%                 fineGridFactor = 10;
%             case 'previous'
%                 fineGridFactor = 10;
            case 'pchip'
                fineGridFactor = 6;
            case 'cubic'
                fineGridFactor = 6;
            case 'v5cubic'
                fineGridFactor = 6;
            case 'spline'
                fineGridFactor = 3;
            otherwise
                fineGridFactor = 10;
        end
    end
end