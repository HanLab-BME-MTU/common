function [block,Bleng]=findBlock(TS,minLength)
%Find blocks of consecutive intergers in TS that have minimum length of "minLength" 
%
%Usage:
%       [block,Bleng]=findBlock(TS,minLength)
%
%Output
%   block - cell array of the blocks
%   Bleng - length of each block
%Marco Vilela, 1-20-11

if nargin < 2
    
    minLength = 1;
    
end

TS = TS(:);

test       = [TS(1)-2;TS;TS(end)+2];
testDiff   = find( diff( diff( test ) == 1 ) );
%Blocks with length >= 2
testDiff   = reshape( testDiff,2,numel(testDiff)/2);
%Block with length = 1
smallBlock = setdiff( find(diff( test ) ~= 1),testDiff(:) );smallBlock([1 end]) = [];
%Temporal indexes 
[~,idx]    = sort([testDiff(1,:) smallBlock']);
%Convert into cell array
smallBlock = num2cell(smallBlock);
bigBlock   = cellfun( @(x) x(1):x(end),num2cell( testDiff,1 ),'Unif',0);
%Set the right temporal chain of events
block      = [bigBlock smallBlock'];
block      = block(idx); 
% Getting rid of blocks with length < minLength
Bleng      = cell2mat(cellfun(@length,block,'Unif',0));
out        = Bleng < minLength;
Bleng(out) = [];
block(out) = [];