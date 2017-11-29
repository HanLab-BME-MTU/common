function [ vq ] = interpcheb1( varargin )
%interpcheb1 Interpolate Chebyshev polynomial, wrapper for interpft1
%
% INPUT / OUTPUT same as interpft1
%
% x defines the bounded range for xq, which is used to map xq to [-pi,pi]
%
% See also interpft1, interp1

    [x,v,xq,method,fineGridFactor,legacy] = parseinputs(varargin{:});
    
    colon = {':'};
    v_colon = colon(ones(ndims(v)-1,1));
    
    % Mirror v in the first dimension
    v = [v ; v(end-1:-1:2,v_colon{:})];
    
    % Map abscissa into [-pi pi]
    xq = acos((xq - x(1))/(x(2)-x(1))*2-1);
    x = [-pi pi];
    vq = interpft1(x,v,xq,method,fineGridFactor,legacy);
    
end

% Copied verbatim from interpft1/parseinputs
% Replace with current interpft1/parseinputs if needed
function [x,v,xq,method,fineGridFactor,legacy] = parseinputs(varargin)
    % Optional arguments
    method = 'pchip';
    fineGridFactor = [];
    x = [];
    legacy = false;
    if(islogical(varargin{end}) && isscalar(varargin{end}))
        legacy = varargin{end};
        varargin = varargin(1:end-1);
    end
    % Decide whether optional arguments are specified, based existence of x
    switch(length(varargin))
        case 6
            [x, v,xq,method,fineGridFactor,legacy] = varargin{:};
        case 5
            % All args specified
            % x v xq method fineGridFactor
            [x, v,xq,method,fineGridFactor] = varargin{:};
        case 4
            if(ischar(varargin{4}))
                % x v xq method
                [x,v,xq,method] = varargin{:};
            else
                % v xq method fineGridFactor
                [v,xq,method,fineGridFactor] = varargin{:};
            end
        case 3
            if(ischar(varargin{3}))
                % v xq method
                [v, xq, method] = varargin{:};
            else
                % x v xq
                [x, v, xq] = varargin{:};
            end
        case 2
            % v xq

            [v, xq]  = varargin{:};
        otherwise
            error('interpft1:nargin','Incorrect number of inputs');
    end
    % Row values are transposed automatically like interp1
    if(isrow(v))
        v = v(:);
    end
    % Query points are transposed automatically like interp1
    if(isrow(xq))
        if(iscolumn(v) || legacy)
            xq = xq(:);
        end
    end
    % Default x specifying periodic boundary f(x) == f(size(v,1)+x)
    switch(numel(x))
        case 2
        case 0
            x = [1 size(v,1)+1];
        otherwise
            error('interpft1:IncorrectX', ...
            'x must either be empty or have 2 elements.');
    end
    if(~legacy)
        % Not in legacy mode
        xq_sz = size(xq);
        v_sz = size(v);
        same_sz = xq_sz(2:end) == v_sz(2:end);
        % Append ones to right side of v_sz so that
        %    length(xq_sz) == length(v_sz)
        v_sz(length(v_sz)+1:length(xq_sz)) = 1;
        % xq should either be a column such that all points are queried against all trig polynomials
        %    outSz = [length(xq) v_sz(2:end)]
        % OR be of the same size as v for every dimension except the first.
        %    outSz = size(xq)
 
        assert(iscolumn(xq) || all(xq_sz(2:end) == v_sz(2:end)) || all(xq_sz([false ~same_sz]) == 1), ...
                'interpft1:InputDimensions', ...
                'xq should either be a column or have similar dimensions as v if not in legacy mode');
    % else
        % Legacy mode
        % outSz = size(xq) % if v is a vector
        % outSz = [length(xq) v_sz(2:end)] % if xq is a vector, and v is not
        % outSz = [xq_sz v_sz(2:end)] % if xq and v both are not vectors
    end
end