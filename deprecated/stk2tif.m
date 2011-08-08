function stk2tif(filename)
% stk2tif converts STK to TIF using stackRead
%
% SYNOPSIS   stk2tif(filename)
%
% INPUT      filename : string containing the common name of the sequence
%            of .tif files
%
% OUTPUT     none     :    
% 
% REMARKS    files are written to disk into user-selected directory
%       
%
%
% DEPENDENCES   stk2tif uses {stackRead}
%               stk2tif is used by {}
%
% Sylvain Berlemont, 20th Jan 2009

fprintf(2,['Warning: ''' mfilename ''' is deprecated and should no longer be used. Use stk2tiffDirs instead.\n']);

if nargin~=1 || isempty(filename)
    error('Please enter a valid (common) file name for the output files');
end

stack = stackRead;

n = size(stack, 3);
L=length(num2str(n)); 
strg=sprintf('%%.%dd', L);

% Select a directory where the output .tif files will be written
path=uigetdir('','Select output directory');
if path==0
    disp('Aborting...');
    return
end    

h=waitbar(0,'Please wait! Writing .tif files');
for z=1:n
    indxStr=sprintf(strg,z);
    imwrite(stack(:, :, z),[path,filesep,filename,indxStr,'.tif']);
    waitbar(z/n,h);
end
close(h);