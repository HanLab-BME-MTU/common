function [vq] = interpft1(varargin)
% INTERPFT1 Evaluate 1D Fourier series at off-grid points, xq
%
% INPUT
% x - (optional) interval that the periodic function is defined on
%     default: [1 size(v,1)+1]
% v - known values at regular intervals corresponding to
%     linspace(x(1),x(2)-x(1),size(v,1))
% xq - abisssa to determine value of Fourier series at
% method - (optional) interpolation method
%          See interp1 for methods, also available:
%          horner - Use Horner's method from chebfun
%          mmt - matrix multiplication method
%          default: pchip
% fineGridFactor - (optional) Number of grid points to expand to using interpft
%                  default: Depends on method
%                  6 for cubic, pchip, v5cubic
%                  3 for spline
%                  10 for everything else
%
% horner and mmt are very accurate methods, but are slow
% interp1 methods are typically faster but less accurate
%
% OUTPUT
% vq - values of Fourier series at xq

% Mark Kittisopikul
% July 20th, 2016

% TODO:
% 1. Extract mapping code, consider wraparound(N)
% 2. Do horner and mmt need wrapping?

    [x,v,xq,method,fineGridFactor] = parseinputs(varargin{:});
    switch(method)
        case 'horner'
            % use external package if availabile
            str = which('trigtech.horner');
            if(isempty(str))
                warning('interpft1:chebfun_missing', ...
                    'trigtech.horner in chebfun is not available. Using mmt method instead');
                % Do mmt, repeated code
                xq = xq - x(1);
                xq = xq./(x(2) - x(1))*2*pi;
                vq = matrixMultiplicationTransform(v,xq);
                return;
            else
                % trigtech.horner is aligned on the domain [0,2)
                % The output is not normalized by the number of points
                xq = xq - x(1);
                xq = xq./(x(2) - x(1))*2;
                vq = trigtech.horner(xq(:),fftshift(fft(v)),isreal(v))./size(v,1);
                v_sz = size(v);
                xq_sz = size(xq);
                if(xq_sz(2) == 1)
                    xq_sz = xq_sz(1);
                end
                vq = reshape(vq,[xq_sz v_sz(2:end)]);
                return;
            end
        case 'mmt'
            xq = xq - x(1);
            xq = xq./(x(2) - x(1))*2*pi;
            vq = matrixMultiplicationTransform(v,xq);
            return;
    end
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
    % Optional arguments
    x = [];
    method = [];
    fineGridFactor = [];
    % Decide whether optional arguments are specified, based existence of x
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
            % All args specified
            % x v xq method fineGridFactor
            [x, v,xq,method,fineGridFactor] = varargin{:};
        otherwise
            error('interpft1:nargin','Incorrect number of inputs');
    end
    % Row values are transposed automatically like interp1
    if(isrow(v))
        v = v.';
    end
    % Query points are transposed automatically like interp1
    if(isrow(xq))
        xq = xq.';
    end
    % Default x specifying periodic boundary f(x) == f(size(v,1)+x)
    if(isempty(x))
        x = [1 size(v,1)+1];
    else
        assert(numel(x) == 2,'interpft1:IncorrectX', ...
            'x must either be empty or have 2 elements.');
    end
    % Use pchip which is reasonably fast and uses a modest amount of memory
    if(isempty(method))
        method = 'pchip';
    end
    % Courser methods should use a finer grid if none is specified
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

function vq = matrixMultiplicationTransform(v,xq)
%matrixMultiplicationTransform
%
% Adapted from interpft_extrema
    s = size(v);
    scale_factor = s(1);

    % Calculate fft and nyquist frequency
    v_h = fft(v);
    nyquist = ceil((s(1)+1)/2);

    % If there is an even number of fourier coefficients, split the nyquist frequency
    if(~rem(s(1),2))
        % even number of coefficients
        % split nyquist frequency
        v_h(nyquist,:) = v_h(nyquist,:)/2;
        v_h = v_h([1:nyquist nyquist nyquist+1:end],:);
        v_h = reshape(v_h,[s(1)+1 s(2:end)]);
    end
    % Wave number, unnormalized by number of points
    freq = [0:nyquist-1 -nyquist+1:1:-1]';
    
    % calculate angles multiplied by wave number
    theta = bsxfun(@times,xq,shiftdim(freq,-ndims(xq)));
    % waves
    waves = exp(1i*theta);


    % evaluate each wave by fourier coeffient
    % theta and waves have one more dimension than xq, representing
    % frequency
    ndims_waves = ndims(waves); % ndims(xq) + 1 ?
    % permute v_h such that it is a 1 by (array dim of fourier series) by
    %                               length of fourier series
    dim_permute = [ndims_waves 2:ndims_waves-1 1];
    
    % sum across waves weighted by Fourier coefficients
    % normalize by the the number of Fourier coefficients
    vq = real(sum(bsxfun(@times,waves,permute(v_h,dim_permute)),ndims_waves))/scale_factor;
end
