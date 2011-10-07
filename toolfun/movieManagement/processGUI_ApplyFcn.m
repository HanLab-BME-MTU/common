function processGUI_ApplyFcn(hObject, eventdata, handles,funParams)
%processGUI_ApplyFcn is a callback called when setting concrete process GUIs
%
%
% Sebastien Besson May 2011 (last modified Oct 2011)

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

% Get the main figure userData
userData = get(handles.figure1, 'UserData');

% Check if the current process is equal to the package process (to cover
% empty processes as well as new subclass processes)
if ~isequal(userData.crtPackage.processes_{userData.procID},userData.crtProc)
    
    if isempty(userData.crtPackage.processes_{userData.procID})
        % Create a new process and set it in the package
        userData.MD.addProcess(userData.crtProc);
        userData.crtPackage.setProcess(userData.procID,userData.crtProc);
    else
        userData.MD.replaceProcess(userData.crtPackage.processes_{userData.procID},userData.crtProc);
    end

       
    % Set font weight of process name bold
    set(userData.handles_main.(['checkbox_' num2str(userData.procID)]),...
            'FontWeight','bold');
end

 % Override the parameters with the GUI set-up ones
parseProcessParams(userData.crtProc,funParams);
 
% ----------------------Sanity Check (II, III check)----------------------

% Do sanity check - only check changed parameters
[status procEx] = userData.crtPackage.sanityCheck(false,'all');

% Draw some bugs on the wall 
for i=find(status)
    userfcn_drawIcon(userData.handles_main,'pass',i,'Current step was processed successfully', true);
end

for i=find(~status)
    userfcn_drawIcon(userData.handles_main,'clear',i,'', true);
end

validProcEx = find(~cellfun(@isempty,procEx));
for i = validProcEx
    % Draw warning label on the i th process
    userfcn_drawIcon(userData.handles_main,'warn',i,...
       sprintf('%s\n',procEx{i}(:).message), true)
end

userData_main = get(userData.mainFig, 'UserData');

if get(handles.checkbox_applytoall, 'Value'),
    moviesId = setdiff(1:numel(userData_main.MD),userData_main.id);
else
    moviesId=[];
end

% Apply setting to all movies
for i = moviesId
    
    % if process classes differ, create a new process with default parameters
    if ~strcmp(class(userData_main.package(i).processes_{userData.procID}),...
            class(userData.crtProc))
        newProcess = userData.procConstr(userData_main.MD(i), ...
            userData_main.package(i).outputDirectory_);
        
        % if package process is empty, add new process and associate it
        if isempty(userData_main.package(i).processes_{userData.procID})
            userData_main.MD(i).addProcess(newProcess);
            userData_main.package(i).setProcess(userData.procID, newProcess);
        else
            userData_main.MD(i).replaceProcess(userData_main.package(i).processes_{userData.procID},newProcess);
        end
    end
    
    % Override the parameters with the GUI defeined
    parseProcessParams(userData_main.package(i).processes_{userData.procID},...
        funParams);
    
    % Do sanity check - only check changed parameters
    [status procEx] = userData_main.package(i).sanityCheck(false,'all');
    
    for j=find(status)
        userData_main.statusM(i).IconType{j} = 'pass';
        userData_main.statusM(i).Msg{j} ='Current step was processed successfully';
    end
    
    for j=find(~status)
        userData_main.statusM(i).IconType{j} = 'clear';
        userData_main.statusM(i).Msg{j} ='Current step was processed successfully';
    end

    % Draw some bugs on the wall
    validProcEx = find(~cellfun(@isempty,procEx));
    for j = validProcEx
        % Record the icon and message to user data
        userData_main.statusM(i).IconType{j} = 'warn';
        userData_main.statusM(i).Msg{j} = sprintf('%s\n',procEx{j}(:).message);
    end
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