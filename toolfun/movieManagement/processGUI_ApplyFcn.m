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

 % Override the parameters with the GUI set-up ones
parseProcessParams(userData.crtProc,funParams);
 
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
       sprintf('%s\n',procEx{i}(:).message), true)
end
% Refresh user data !!
userData_main = get(userData.mainFig, 'UserData');

% -------------------- Apply setting to all movies ------------------------
if get(handles.checkbox_applytoall, 'Value')

    moviesId = setdiff(1:numel(userData_main.MD),userData_main.id);
    for i = moviesId
                
        % if package process is empty, create a new process with default 
        % parameters and add to MovieData and package's process list
        if isempty(userData_main.package(i).processes_{userData.procID})            
            newProcess = userData.procConstr(userData_main.MD(i), ...
                userData_main.package(i).outputDirectory_);
            userData_main.MD(i).addProcess(newProcess);
            userData_main.package(i).setProcess(userData.procID, newProcess);
        end

        % Override the parameters with the GUI defeined
        parseProcessParams(userData_main.package(i).processes_{userData.procID},...
            funParams);
        
        % Do sanity check - only check changed parameters
        procEx = userData_main.package(i).sanityCheck(false,'all');
        
        % Draw some bugs on the wall
        validProcEx = find(~cellfun(@isempty,procEx));
        for j = validProcEx
            % Record the icon and message to user data
            userData_main.statusM(i).IconType{j} = 'warn';
            userData_main.statusM(i).Msg{j} = sprintf('%s\n',procEx{j}(:).message);
        end
    end
    
    % Save user data
    set(userData.mainFig, 'UserData', userData_main)
end
% -------------------------------------------------------------------------

% Store the applytoall choice for this particular process
userData_main.applytoall(userData.procID)=get(handles.checkbox_applytoall,'Value');
% Save user data
set(userData.mainFig, 'UserData', userData_main)
set(handles.figure1, 'UserData', userData);
guidata(hObject,handles);
delete(handles.figure1);
end