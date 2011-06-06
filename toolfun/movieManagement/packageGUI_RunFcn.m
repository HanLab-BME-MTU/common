function packageGUI_RunFcn(hObject,eventdata,handles)
% This is a common section of code called by pushbutton_run_Callback
% when user click the "Run" button on package control panels.
%
% Chuangang Ren 11/2010
% Sebastien Besson 5/2011

ip = inputParser;
ip.addRequired('hObject',@ishandle);
ip.addRequired('eventdata',@(x) isstruct(x) || isempty(x));
ip.addRequired('handles',@isstruct);
ip.parse(hObject,eventdata,handles);

userData = get(handles.figure1,'UserData');

movieRun = []; % id of movie to run
numMD = length(userData.MD); % number of movies


% Get check box status of current movie and update user data
userData.statusM(userData.id).Checked = userfcn_saveCheckbox(handles);
set(handles.figure1, 'UserData', userData)

% Determine the movie(s) to be processed
if ~get(handles.checkbox_runall, 'Value')
    
    if any(userData.statusM(userData.id).Checked)
        movieRun = cat(2, movieRun, userData.id);
    end
else
    
    % Check other movies
    for i = 0:numMD-1
        k = mod( userData.id + i, numMD);
        if k == 0
            k = numMD;
        end
        
        if any(userData.statusM(k).Checked)
            movieRun = cat(2, movieRun, k);
        end
    end
end

if isempty(movieRun)
    warndlg('No step is selected, please select a step to process.','No Step Selected','modal');
    return
end

% ----------------------- Pre-processing examination ----------------------


% movie exception (same length of movie data)
movieException = cell(1, numMD);

procCheck = cell(1, numMD);%  id of checked processes
procRun = cell(1, numMD);%  id of processes to run
optionalProcID = cell(1, numMD);% id of first-time-run optional process

for x = movieRun
    
    
    procCheck{x} = find(userData.statusM(x).Checked);
    
    % Check if process exist
    for i = procCheck{x}
        if isempty (userData.package(x).processes_{i})
            
            ME = MException('lccb:run:setup', 'Step %d is not set up yet. Tip: when step is set up successfully, the step name becomes bold.',i);
            movieException{x} = cat(2, movieException{x}, ME);
            
        end
    end
    
    if ~isempty(  movieException{x} )
        continue
    end
    
    % Check if selected processes have alrady be successfully run
    % If force run, re-run every process that is checked
    if ~get(handles.checkbox_forcerun, 'Value')
        
        k = true;
        for i = procCheck{x}
            
            if  ~( userData.package(x).processes_{i}.success_ && ...
                    ~userData.package(x).processes_{i}.procChanged_ ) || ...
                    ~userData.package(x).processes_{i}.updated_
                
                k = false;
                procRun{x} = cat(2, procRun{x}, i);
            end
        end
        if k
            movieRun = setdiff(movieRun, x);
            continue
        end
    else
        procRun{x} = procCheck{x};
    end
    
    
    
    % Package full sanity check. Sanitycheck every checked process
    procEx = userData.package(x).sanityCheck(true, procRun{x});
    
    % Return user data !!!
    set(handles.figure1, 'UserData', userData)
    
    for i = procRun{x}
        if ~isempty(procEx{i})
            
            % Check if there is fatal error in exception array
            if strcmp(procEx{i}(1).identifier, 'lccb:set:fatal') || ...
                    strcmp(procEx{i}(1).identifier, 'lccb:input:fatal')
                
                % Sanity check error - switch GUI to the x th movie
                if x ~= userData.id
                    set(handles.popupmenu_movie, 'Value', x)
                    % Quick fix for callback incompatibility betwen packageGUI and
                    % oldpackageGUIs - to be solved before release
                    stack = dbstack;
                    if strcmp(stack(4).name,'packageGUI')
                        switchMovie_Callback(handles.popupmenu_movie, [], handles) % user data retrieved, updated and submitted
                    else
                        popupmenu_movie_Callback(handles.popupmenu_movie, [], handles) % user data retrieved, updated and submitted
                    end
                end
                
                userfcn_drawIcon(handles,'error', i, procEx{i}(1).message, true); % user data is retrieved, updated and submitted
                
                ME = MException('lccb:run:sanitycheck', 'Step %d %s: \n%s', i,userData.package(x).processes_{i}.name_, procEx{i}(1).message);
                movieException{x} = cat(2, movieException{x}, ME);
                
            end
        end
    end
    
    % Refresh user data !!!
    userData = get(handles.figure1, 'UserData');
    
    
end

% --------------------- pre-processing examination ends -------------------

% Ok, now all evils are in movieException (1 x movielength  cell array), if there is any
% if yes - abort program and popup a error report
% if no - continue to process movie data
if isempty(movieRun)
    warndlg('All selected steps have been processed successfully. Please check the ''Force Run'' check box if you want to re-process the successful steps.','No Step Selected','modal');
    return
end

temp = find(~cellfun(@(x)isempty(x), movieException, 'UniformOutput', true));

if ~isempty(temp)
    msg = [];
    for i = 1:length(temp)
        if i == 1
            msg = strcat(msg, sprintf('Movie %d - %s:', temp(i), userData.MD(temp(i)).movieDataFileName_));
        else
            msg = strcat(msg, sprintf('\n\n\nMovie %d - %s:', temp(i), userData.MD(temp(i)).movieDataFileName_));
        end
        for j = 1:length(movieException{temp(i)})
            msg = strcat(msg, sprintf('\n-- %s', movieException{temp(i)}(j).message));
        end
        
    end
    msg = strcat(msg, sprintf('\n\n\nPlease solve the above problems before continuing. The Movie(s) couldn’t be processed.'));
    titlemsg = sprintf('Processing could not be continued for the following reasons:');
    
    % if msgboxGUI exist
    if isfield(userData, 'msgboxGUI') && ishandle(userData.msgboxGUI)
        delete(userData.msgboxGUI)
    end
    
    userData.msgboxGUI = msgboxGUI('title',titlemsg,'text', msg);
    return
    
end

% ------------------------ Start Processing -------------------------------
kk = 0;
for x = movieRun
    
    kk = kk+1;
    if x ~= userData.id
        
        set(handles.popupmenu_movie, 'Value', x)
        set(handles.figure1, 'UserData', userData)
        % Update the movie pop-up menu in the main package GUI
        packageGUI('switchMovie_Callback',handles.popupmenu_movie, [], handles)
        userData = get(handles.figure1, 'UserData');
    end
    
    % Find first-time-run optional process ID
    for i = intersect(procRun{x}, userData.optProcID);
        if ~userData.package(x).processes_{i}.success_
            optionalProcID{x} = cat(2, optionalProcID{x}, i);
        end
    end
    
    % Clear icons of selected processes
    % Return user data !!!
    set(handles.figure1, 'UserData', userData)
    userfcn_drawIcon(handles,'clear',procRun{x},'',true); % user data is retrieved, updated and submitted
    % Refresh user data !!!
    userData = get(handles.figure1, 'UserData');
    
    % Disable 'Run' button
    set(handles.pushbutton_run, 'Enable', 'off')
    set(handles.checkbox_forcerun, 'Enable', 'off')
    set(handles.checkbox_runall, 'Enable', 'off')
    set(handles.text_status, 'Visible', 'on')
    
    % Run algorithms!
    try
        % Return user data !!!
        set(handles.figure1, 'UserData', userData)
        
        for i = procRun{x}
            set(handles.text_status, 'String', sprintf('Step %d - Processing %d of %d movies total ...', i, kk, length(movieRun)) )
            userfcn_runProc_dfs(i, procRun{x}, handles); % user data is retrieved, updated and submitted
            
        end
        
    catch ME
        
        % Save the error into movie Exception cell array
        movieException{x} = ME;
        
        procRun{x} = procRun{x}(procRun{x} < i);
        optionalProcID{x} = optionalProcID{x}(optionalProcID{x} < i); 
    end
    
    % Refresh user data !!!
    userData = get(handles.figure1, 'UserData');
    set(handles.pushbutton_run, 'Enable', 'on')
    set(handles.checkbox_forcerun, 'Enable', 'on')
    set(handles.checkbox_runall, 'Enable', 'on')
    set(handles.text_status, 'Visible', 'off')
    
    % ------- Check optional processes ----------
    
    % Return user data !!!
    set(handles.figure1, 'UserData', userData)
    % In here, optionalProcID are successfuly first-time-run optional process ID
    if ~isempty(optionalProcID{x})
        
        procEx = userData.crtPackage.checkOptionalProcess(procRun{x}, optionalProcID{x});
        
        for i = 1:size(userData.dependM, 1)
            if ~isempty(procEx{i})
                
                userfcn_drawIcon(handles,'warn',i,procEx{i}(1).message, true); % user data is retrieved, updated and submitted
                
            end
        end
    end
    
end

% ----------------------------- Create error report ---------------------------------------

temp = find(~cellfun(@(x)isempty(x), movieException, 'UniformOutput', true));

if ~isempty(temp)
    msg = [];
    for i = 1:length(temp)
        if i == 1
            msg = strcat(msg, sprintf('Movie %d - %s:\n\n%s', ...
                temp(i), userData.MD(temp(i)).movieDataFileName_, movieException{x}.message));
        else
            msg = strcat(msg, sprintf('\n\n\nMovie %d - %s:\n\n%s', ...
                temp(i), userData.MD(temp(i)).movieDataFileName_, movieException{x}.message));
        end
    end
    
    msg = strcat(msg, sprintf('\n\n\nPlease verify your settings are correct. Feel free to contact us if you have question regarding this error.\n\nPlease help us improve the software by clearly reporting the scenario when this error occurs, and the above error information to us (error information is also displayed in Matlab command line).\nFor contact information please refer to the following URL:\n\nlccb.hms.harvard.edu/software.html'));
    
    % if msgboxGUI exist
    if isfield(userData, 'msgboxGUI') && ishandle(userData.msgboxGUI)
        delete(userData.msgboxGUI)
    end
    if length(temp) == 1
        userData.msgboxGUI = msgboxGUI('title','The processing of following movie is terminated by run time error:','text', msg);
    else
        userData.msgboxGUI = msgboxGUI('title','The processing of following movies are terminated by run time errors:','text', msg);
    end
    
else
    if length(movieRun) > 1
        successMsg = 'All your movies have been processed successfully.';
    else
        successMsg = 'Your movie has been processed successfully.';
    end
    
    userData.iconHelpFig = helpdlg(successMsg, [userData.crtPackage.name_ 'Package']);
    set(handles.figure1, 'UserData', userData)
end
end

function userfcn_runProc_dfs (i, procRun, handles)  % throws exception

% Set user Data
userData = get(handles.figure1, 'UserData');

parentRun = [];
parentIndex = find(userData.crtPackage.depMatrix_(i,:));

% if current process i have dependency processes    
if ~isempty(parentIndex)  
    for j = parentIndex
        % if parent process is one of the processes need to be run
        % if parent process has already run successfully
        if any(j == procRun) && ~userData.crtPackage.processes_{j}.success_
            parentRun = horzcat(parentRun,j); %#ok<AGROW>
        end
    end
    % if above assumptions are yes, recursively run parent process' dfs fcn
    if ~isempty(parentRun)
        for j = parentRun
            userfcn_runProc_dfs (j, procRun, handles)
        end
    end
end

try
    userData.crtPackage.processes_{i}.run(); % throws exception
catch ME
    rethrow(ME) %%%%
    % Determine if it is an unexpected error
%     idSplit = regexp(ME.identifier, ':', 'split');
%     
%     if isempty(idSplit{1}) || strcmp(idSplit{1}, 'lccb')
%         errorText = sprintf('Step %d - %s: Runtime error \n%s',i, userData.crtPackage.processes_{i}.name_, ME.message);
%     else
%     
%         errorText = sprintf...
%         ('Step %d - %s: Unexpected runtime error \nIdentifier: %s\nMessage: %s\nErrorfcn: %s\nErrorline: %u',...
%         i, userData.crtPackage.processes_{i}.name_, ME.identifier, ME.message, ME.stack(1).name, ME.stack(1).line);
%         
%         display(sprintf('\n??? %s', errorText))
%     end
%     
%     
% 
%     set(handles.figure1, 'UserData', userData)
%     userfcn_drawIcon(handles,'error',i,errorText, true); % user data is retrieved, updated and submitted
%     userData = get(handles.figure1, 'UserData');
%     
%     ME2 = MException('lccb:runtime:fatal', errorText);
%     ME2 = addCause(ME2, ME);
%     throw(ME2);
end

% After successfully processed, determine if dependent processes are updated.
% If no, set current process updated = false, and draw warning icon
% if yes, set current process updated = true, and draw pass icon
l = true;

% Return user data !!!
set(handles.figure1, 'UserData', userData)

for k = parentIndex
   if ~userData.crtPackage.processes_{k}.updated_ 
           
       userData.crtPackage.processes_{i}.setUpdated(false);
       userfcn_drawIcon(handles,'warn',i,...
         ['Current step is processed successfully. But it is found to be out of date.'...
              'Please make sure the dependent steps are up to date.'], true); % user data is retrieved, updated and submitted
      l = false;
       break
   end
end

if l
    userData.crtPackage.processes_{i}.setUpdated(true);
    userfcn_drawIcon(handles,'pass',i,...
                                'Current step is processed successfully', true); % user data is retrieved, updated and submitted
end

set(handles.(['pushbutton_show_',num2str(i)]),'Enable','on');
end