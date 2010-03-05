function fixImageNumbers(directory)
% 
% fixImageNumbers(directory)
% 
% This function goes through every numbered image file in the specified
% directory and renames files which have an ending like _t1, _t2 etc. To
% _t01, _t02 so that dir/ls returns them in the correct order.
% 
% Hunter Elliott, 11/2009
%

stackFiles = imDir(directory);

nStack = length(stackFiles);

nDig = floor(log10(nStack))+1;
fString = ['%0' num2str(nDig) '.0f'];


if nStack > 0
    
    disp(['Renaming ' num2str(nStack) ' image files...'])
    
    for i = 1:nStack
        
        %Get the index of the last t
        iT = max(regexp(stackFiles(i).name,'t'));
        
        %And of the file extension
        iLFS = max(regexp(stackFiles(i).name,'\.'));
        
        iFrame = str2double(stackFiles(i).name(iT+1:iLFS-1));
        
        newName = [stackFiles(i).name(1:iT) num2str(iFrame,fString) stackFiles(i).name(iLFS:end)];
        
        if ~strcmp(newName,stackFiles(i).name)
            movefile([directory filesep stackFiles(i).name],...
                [directory filesep newName]);                
        end
    end
else 
    disp('No image files found in specified directory...')
    return
end

disp('Finished!')




