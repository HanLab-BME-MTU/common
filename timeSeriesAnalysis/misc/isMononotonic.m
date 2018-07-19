function out = isMononotonic(TS,varargin)
% This function outputs 1 if TS is monotonic and 0 otherwise
%
% USAGE: out = isMononotonic(TS,'type',0,'direction',1)
%
% INPUT:
%       TS        : vector of points
%       type      : scalar - 0 for monotonic
%                            1 for strictly monotonic
%       
%       direction :  1 for increasing
%                   -1 for decrasing
%
% OUTPUT:
%       out: 1 if TS satisfy the test and 0 otherwise
%
%
% Marco Vilela, 2012


%% Parsing the input ******************************************************
ip = inputParser;
ip.addRequired('TS',@(x) isvector(x));
ip.addOptional('type',0,@isscalar);
ip.addOptional('direction',1,@isscalar);

ip.parse(TS,varargin{:});
type      = ip.Results.type;
direction = ip.Results.direction;

%% 
if type 
    
    test = @gt;
    
else
    
    test = @ge;
    
end

out = all( test( diff( direction*TS ), 0 ) );