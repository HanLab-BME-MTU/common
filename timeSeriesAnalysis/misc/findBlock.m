function [block,Bleng]=findBlock(S,min_len)
%[block,Bleng]=find_block(S)
%Find blocks of consecutive intergers in S that have minimum length of "min_len" 
%Output
%   block - cell array of the blocks
%   Bleng - length of each block
%Marco Vilela, 1-20-11
if nargin < 2
    min_len=1;
end

n=length(S);
cc1=1;
bk=1;
block=[];
Bleng=0;
while cc1 < n
    cc=1;
    flag=1;
    while flag
        if (diff(S(cc1+cc-1:cc1+cc)) ~= 1) || (cc1+cc) == n
            flag=0;
            break;
        end
        cc=cc+1;
    end
    if cc >= min_len
        block{bk}=S(cc1:cc1+cc-1);
        Bleng(bk)=cc;
        bk=bk+1;
    end
    cc1=cc1+cc;
end

if ~isempty(block)
    if (S(end) - block{bk-1}(end)) == 1
        block{bk-1}(cc+1)=S(end);
        Bleng(bk-1)=Bleng(bk-1)+1;
        for i=1:bk-1
            block{i}=reshape(block{i},length(block{i}),1);
        end

    elseif min_len == 1
        block{bk}=S(end);
        Bleng(bk)=1;
    end
elseif n == 1
     block{bk}=S(end);
     Bleng(bk)=1;
end

%Make sure it always return a column vector
nB = length(block);
for i =1:nB
    [nP,nV] = size(block{i});
    if nV > nP
        block{i} = block{i}';
    end
end

