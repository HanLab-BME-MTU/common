function packageGUI_RefreshFcn(handles, type)
% GUI tool function: this function is called by movie explorer when 
% switching between differenct movies. 
% 
% Input: 
%       handles - the "handles" of package GUI control panel
%       type - 'initialize': used when movie is loaded to GUI for the first time
%              'refresh': used when movie had already been loaded to GUI
%
%
% Chuangang Ren 08/2010

% Input check
ip = inputParser;
ip.addRequired('handles',@isstruct);
ip.addRequired('type',@(x) any(strcmp(x,{'initialize','refresh'})));
ip.parse(handles,type)


userData = get(handles.figure1, 'UserData');
l = size(userData.dependM, 1);
k = zeros(1,l);

% Reset GUI
userfcn_drawIcon(handles, 'clear', 1:l);
userfcn_enable(1:l, 'on', handles)
for i = 1:l
    set(handles.(['checkbox_',num2str(i)]),'FontWeight','normal','Value',0);
    set(handles.(['pushbutton_show_',num2str(i)]),'Enable','off');
end

% Set movie data path
set(handles.edit_path, 'String', ...
    [userData.MD(userData.id).movieDataPath_ filesep userData.MD(userData.id).movieDataFileName_ ])


% ----------------------------- Initialize --------------------------------

if strcmp(type, 'initialize')

    % Package Sanity Check
    [status procEx] = userData.crtPackage.sanityCheck(true, 'all');

    % Draw successful processes
    for i=find(status)
        userfcn_drawIcon(handles,'pass',i,'Current step was processed successfully', true);
    end
    
    % Clear unsuccesful processes
    for i=find(~status)
        userfcn_drawIcon(handles,'clear',i,'', true);
    end
    
    % Draw warnings
    validProcEx = find(~cellfun(@isempty,procEx));
    for i = validProcEx
        if strcmp(procEx{i}(1).identifier, 'lccb:set:fatal')
            statusType='error';
        else
            statusType='warn';
        end
        userfcn_drawIcon(handles,statusType,i,...
            sprintf('%s\n',procEx{i}(:).message), true);
    end


% ----------------------------- Refresh -----------------------------------    
elseif strcmp(type, 'refresh')
    for i = 1: l
    
       % Draw icons
       if ~isempty(userData.statusM(userData.id).IconType{i})
            userfcn_drawIcon(handles, userData.statusM(userData.id).IconType{i}, i, userData.statusM(userData.id).Msg{i}, false);
       end

    end
end  
% -------------------------------------------------------------------------

for i = 1: l
    
    % If process is checked, check and enable the process and enable decendent
    % processes
    if userData.statusM(userData.id).Checked(i)
        k(i) = 1;
        set(handles.(['checkbox_',num2str(i)]),'Value',1,'Enable','on');
        userfcn_lampSwitch(i, 1, handles)
    end
    
    % Bold the Name of Existing Process
    
    if ~isempty(userData.crtPackage.processes_{i})
        set(handles.(['checkbox_',num2str(i)]),'FontWeight','bold');
        % Set Up Uicontrols Enable/Disable
        
        % If process's sucess = 1, allow output visualizatoin
        if userData.crtPackage.processes_{i}.success_
            k(i) = 1;
            set(handles.(['pushbutton_show_',num2str(i)]),'Enable','on');
        end
    end
    
end

tempDependM = userData.dependM;
tempDependM(:,logical(k)) = zeros(l, nnz(k));

% Checkbox enable/disable set up
userfcn_enable(find (any(tempDependM==1,2)), 'off',handles);


if strcmp(type, 'initialize')
    userData.statusM(userData.id).Visited = true;
    set(handles.figure1, 'UserData', userData)
end
