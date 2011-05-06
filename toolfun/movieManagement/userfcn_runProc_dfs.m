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