function movieArray = relocateMovieData(parentDir,confirmChange)

% movieArray = relocateMovieData(parentDir)
% 
% This function alters the directory and file locations stored in a set of
% movieDatas that have been moved from one directory to another to reflect
% their new location. Use this function AFTER moving the movies and
% movieDatas to their new location.
% 
% Input: 
% 
% parentDir - The new parent directory containing all the movies &
%                   movieDatas which have been moved. Optional. If not
%                   input, the user will be asked to specify a directory.
% 
% confirmChange - If true, the user will be asked to confirm the directory
%                 change before continuing. If false, the changes will be
%                 made without confirmation. Default is true. 
%                 NOTE: If this option is false, the parent directory MUST
%                 be specified.
% 
% 
% Output:
% 
% movieArray - cell array of the new, relocated movieDatas.
%
% Hunter Elliott, 10/2009
%


%%  ----- Init ------%%


showOutput = true; %Determines whether progress is shown at the command prompt
yesAllStr = 'Yes, same for all movies';

if nargin < 2 || isempty(confirmChange)
    confirmChange = true;
end

if nargin < 1 || isempty(parentDir)    
    if ~confirmChange
        parentDir = uigetdir(pwd,'Please specify the parent directory containing all the movies:');        
    else
       error('If the confirm option is set to false, the parent directory MUST be input!!') 
    end
end




movieArray = [];


%Check the directory
if ~exist(parentDir,'dir')    
    error('Specified parent directory does not exist!')
end


%% --- Find MovieDatas ---%%


if showOutput
    disp('Searching all sub-directories for files named movieData.mat... Please be patient!')
end

%Search all sub-directories for files named movieData.mat
mdFileList = searchFiles('movieData.mat',[],parentDir,1,'all',1);

nMov = length(mdFileList);

if nMov < 1
    error('No movieData(s) found in specified parent directory!')
elseif showOutput
   disp(['Found ' num2str(nMov) ' movieDatas for relocation.'])    
end


%% ---- Relocate MovieDatas ---- %%
%Go through each movieData and compare it's current location to the
%location specified in the movieData, then change the movieData to match

movieArray = cell(nMov,1);
bPressed = 'No';

for iMov = 1:nMov
    
   
    if showOutput
        disp(['Relocating movieData ' num2str(iMov) '...'])
    end

    %Load the current movieData
    movieArray{iMov} = load(mdFileList{iMov},'movieData');
    movieArray{iMov} = movieArray{iMov}.movieData;

    %Remove the actual file name from the specifier to ease comparison
    mdFileList{iMov} = regexprep(mdFileList{iMov},[filesep 'movieData.mat'],'');

    %Convert file seperators in case there has been a change of OS
    movieArray{iMov} = rReplace(movieArray{iMov},'/|\',filesep);        

    if ~strcmp(bPressed,yesAllStr)

        %Compare the old and new locations of the analysis directory to
        %determine the change

        newDir = mdFileList{iMov};
        oldDir = movieArray{iMov}.analysisDirectory;
        
        %Remove any trailing file seperators
        if strcmp(newDir(end),filesep)
            newDir = newDir(1:end-1);
        end
        if strcmp(oldDir(end),filesep)
            oldDir = oldDir(1:end-1);
        end
        
        %Check for the case that the old and new directories are identical
        if strcmp(oldDir,newDir)
            
            %Don't ask - it's identical
            bPressed = 'No';
            disp('Old and new movieData''s are identical!')
            
        else

            %Find the first character where the two directories differ:
            stillMatch = true;
            j=1;
            while stillMatch        
                stillMatch = strncmp(newDir(end:-1:1),oldDir(end:-1:1),j);
                j=j+1;        
            end
            j = j-2; %Remove offset in j.



            % 'Snap' the location of the difference to the next highest directory
            % to correct for directory names which are partially equivalent

            %Find the first seperator they have in common
            iFs = regexp(newDir,filesep);
            iFs = (length(newDir)+1) - iFs(find(iFs > (length(newDir)-j),1));

            %Get the old and new parent directories
            oldPdir = oldDir(1:end-iFs);
            newPdir = newDir(1:end-iFs);

            disp('Old parent directory:')
            disp(oldPdir)
            disp('New parent directory:')
            disp(newPdir)       

            if confirmChange
                %Force the user to confirm this change is correct
                bPressed = questdlg('Are the old and new dirctories correct? (see command prompt)','Confirm Directory Change','Yes',yesAllStr,'No','No');
            else
                bPressed = yesAllStr;
            end
        end
    end
    
    if strcmp(bPressed,'Yes') || strcmp(bPressed,yesAllStr)
    
        movieArray{iMov} = rReplace(movieArray{iMov},regexptranslate('escape',oldPdir),regexptranslate('escape',newPdir));
        updateMovieData(movieArray{iMov})        
        if showOutput
            disp(['Movie ' num2str(iMov) ' successfully relocated!'])
        end
    elseif showOutput
        disp(['Movie ' num2str(iMov) ' NOT relocated - movieData unchanged'])
        if isempty(bPressed)
            isempty(bPressed)
            return
        end
    end
    
end

disp(['Done relocating ' num2str(nMov) ' movies.'])


