function processGUI_ApplyFcn(hObject, eventdata, handles,funParams)
% Common saving of concrete process GUIs settings
%
%
% Sebastien Besson May 2011

% Check input
ip = inputParser;
ip.addRequired('hObject',@ishandle);
ip.addRequired('eventdata',@(x) isstruct(x) || isempty(x));
ip.addRequired('handles',@isstruct);
ip.addRequired('funParams',@isstruct);
ip.parse(hObject,eventdata,handles,funParams);


if get(handles.checkbox_applytoall, 'Value')
    confirmApplytoAll = questdlg(...
        ['You are about to copy the current process settings to all movies.'...
        ' Previous settings will be lost. Do you want to continue?'],...
        'Apply settings to all movies','Yes','No','Yes'); 
    if ~strcmp(confirmApplytoAll,'Yes'),
        set(handles.checkbox_applytoall,'Value',0);
        return
    end
end

% Call back function of 'Apply' button
userData = get(handles.figure1, 'UserData');
userData.crtProc.setPara(funParams);
 
% --------------------------------------------------
userData_main = get(userData.mainFig, 'UserData');
% If this is a brand new process, attach current process to MovieData and 
% package's process list 
if isempty(userData.crtPackage.processes_{userData.procID})
    
    % Add new process to both process lists of MovieData and current package
    userData_main.MD(userData_main.id).addProcess(userData.crtProc);
    userData.crtPackage.setProcess(userData.procID,userData.crtProc);
    
    % Set font weight of process name bold
    set(userData.handles_main.(['checkbox_' num2str(userData.procID)]),...
            'FontWeight','bold');
end

% ----------------------Sanity Check (II, III check)----------------------

% Do sanity check - only check changed parameters
procEx = userData.crtPackage.sanityCheck(false,'all');

% Return user data !!!
set(userData.mainFig, 'UserData', userData_main)

% Draw some bugs on the wall 
validProcEx = find(~cellfun(@isempty,procEx));
for i = validProcEx
    % Draw warning label on the i th process
    userfcn_drawIcon(userData.handles_main,'warn',i,...
        procEx{i}(1).message, true)
end
% Refresh user data !!
userData_main = get(userData.mainFig, 'UserData');

% -------------------- Apply setting to all movies ------------------------
if get(handles.checkbox_applytoall, 'Value')

    for x = 1: length(userData_main.MD)
        
        if x == userData_main.id
            continue
        end
        
        % Customize funParams to other movies
        % OutputDirectory - package output directory      
        funParams.OutputDirectory  = [userData_main.package(x).outputDirectory_  filesep 'masks'];
        % if new process, create a new process with funParas and add to
        % MovieData and package's process list
        if isempty(userData_main.package(x).processes_{userData.procID})            
            process = userData.procConstr(userData_main.MD(x), ...
                userData_main.package(x).outputDirectory_, funParams);
            userData_main.MD(x).addProcess(process);
            userData_main.package(x).setProcess(userData.procID, process);            
        else
            % if process exist, replace the funParams with the new one
            userData_main.package(x).processes_{userData.procID}.setPara(funParams)
        end
        
        
        % Do sanity check - only check changed parameters
        procEx = userData_main.package(x).sanityCheck(false,'all');
        
        % Draw some bugs on the wall
        validProcEx = find(~cellfun(@isempty,procEx));
        for i = validProcEx
            % Record the icon and message to user data
            userData_main.statusM(x).IconType{i} = 'warn';
            userData_main.statusM(x).Msg{i} = procEx{i}(1).message;
        end
    end
    
    % Save user data
    set(userData.mainFig, 'UserData', userData_main)
end
% -------------------------------------------------------------------------

% Save user data
userData_main.applytoall(userData.procID)=get(handles.checkbox_applytoall,'Value');
set(userData.mainFig, 'UserData', userData_main)
set(handles.figure1, 'UserData', userData);
guidata(hObject,handles);
delete(handles.figure1);
end