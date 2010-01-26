function fieldOut = dirFieldConvert(fieldIn)

%Checks if the input field is a string and if it is it converts it
%recursively
% For use with structfun on big structures containing directory names

%Hunter Elliott 2/2009

if ischar(fieldIn) 
    
    if any(regexp(fieldIn,'/')) %If unix
        
        currDrive = getDriveName(fieldIn);
        if ~isempty(currDrive)
            fieldOut = dirUnix2PC(fieldIn,convertDriveName(currDrive));
        else
            fieldOut = fieldIn;
        end
        
    elseif any(regexp(fieldIn,'\')) %If windows
        
        currDrive = getDriveName(fieldIn);
        if ~isempty(currDrive)
            fieldOut = dirPC2Unix(fieldIn,convertDriveName(currDrive));
        else
            fieldOut = fieldIn;
        end
        
        
    else
        fieldOut = fieldIn;
    end
    
    
    
    
elseif isstruct(fieldIn)
    %Otherwise, recursively convert
    fNames = fieldnames(fieldIn);    
    for j = 1:length(fNames)%Call this function on each sub-field
                                
        %Initialize the output
        fieldOut.(fNames{j}) = fieldIn.(fNames{j});
        
        if isstruct(fieldIn.(fNames{j})) %If its a structure, convert each element seperately
            for k = 1:length(fieldIn.(fNames{j}));        
                fieldOut.(fNames{j})(k) = dirFieldConvert(fieldIn.(fNames{j})(k));
            end                        
        else
            fieldOut.(fNames{j}) = dirFieldConvert(fieldIn.(fNames{j}));            
        end
        
    end
        
        
    
else %If it's not a string or a structure, just return it
    
    fieldOut = fieldIn;
    
end