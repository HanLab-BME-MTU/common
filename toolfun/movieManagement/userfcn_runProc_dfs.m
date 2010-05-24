function userfcn_runProc_dfs (i, procRun, MD, handles)  % throws exception

% Set user Data
userData = get(handles.figure1, 'UserData');

parentRun = [];
parentIndex = find(userData.crtPackage.depMatrix_(i,:));

% if current process i have dependency processes    
if ~isempty(parentIndex)  
    for j = parentIndex
        % if parent process is one of the processes need to be run
        % if parent process has already run
        if any(j == procRun) && ~userData.crtPackage.processes_{j}.success_
            parentRun = horzcat(parentRun,j); %#ok<AGROW>
        end
    end
    % if above assumptions are yes, recursively run parent process' dfs fcn
    if ~isempty(parentRun)
        for j = parentRun
            userfcn_runProc_dfs (j,procRun,MD, handles)
        end
    end
end
try
    userData.crtPackage.processes_{i}.runProcess; % throws exception
catch ME
    
    errorText = sprintf...
        ('Runtime error! You may report the following error information to us:\n\n Identifier: %s\n Message: %s\n Errorfcn: %s\n Errorline: %u',...
        ME.identifier,ME.message,ME.stack(1).name,ME.stack(1).line');
    
    userData.crtPackage.processes_{i}.setSuccess(false);
    userfcn_drawIcon(handles,'error',i,errorText);
    
    ME2 = MException('lccb:runtime:fatal', errorText);
        
    ME2 = addCause(ME2, ME);
    throw(ME2);
end

userData.crtPackage.processes_{i}.setSuccess(true);
userData.crtPackage.processes_{i}.setProcChanged(false);

% After successfully processed, determine if dependent processes are updated.
% If no, set current process updated = false, and draw warning icon
% if yes, set current process updated = true, and draw pass icon
l = true;

    for k = parentIndex
       if ~userData.crtPackage.processes_{k}.updated_ 
           userData.crtPackage.processes_{i}.setUpdated(false);
           userfcn_drawIcon(handles,'warn',i,...
             ['Current step is processed successfully. But it is found to be out of date.'...
              'Please make sure the predecessor steps are up to date']);
          l = false;
           break
       end
    end

if l
    userData.crtPackage.processes_{i}.setUpdated(true);
    userfcn_drawIcon(handles,'pass',i,...
                                'Current step is processed successfully');
end

eval([ 'set(handles.pushbutton_show',num2str(i),', ''enable'', ''on'');']);


