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
nProc = size(userData.dependM, 1);

% Reset GUI
% userfcn_drawIcon(handles, 'clear', 1:nProc);
% userfcn_enable(1:nProc, 'on', handles)


% Set movie data path
set(handles.edit_path, 'String', ...
    [userData.MD(userData.id).getPath filesep userData.MD(userData.id).getFilename])

% Bold the Name of Existing Process
setupProc = ~cellfun(@isempty,userData.crtPackage.processes_);
for i=find(setupProc)
    set(handles.(['checkbox_',num2str(i)]),'FontWeight','bold');
end

for i=find(~setupProc)
    set(handles.(['checkbox_',num2str(i)]),'FontWeight','normal');
    set(handles.(['pushbutton_show_',num2str(i)]),'Enable','off');
end


% Run sanityCheck on package 
% if strcmp(type, 'initialize'), full=true; else full=false; end
[status procEx] = userData.crtPackage.sanityCheck(true, 'all');

% Draw successful processes
for i=find(status)
    userfcn_drawIcon(handles,'pass',i,'Current step was processed successfully', true);
    set(handles.(['pushbutton_show_',num2str(i)]),'Enable','on');
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


% If process is checked, check and enable the process and enable decendent
% processes
for i = find(~userData.statusM(userData.id).Checked)
    set(handles.(['checkbox_',num2str(i)]),'Value',0);
end

for i = find(userData.statusM(userData.id).Checked)
     set(handles.(['checkbox_',num2str(i)]),'Value',1,'Enable','on');
     userfcn_lampSwitch(i, 1, handles)
end

% Checkbox enable/disable set up
k=status | userData.statusM(userData.id).Checked;
tempDependM = userData.dependM;
tempDependM(:,logical(k)) = zeros(nProc, nnz(k));
userfcn_enable(find (any(tempDependM==1,2)), 'off',handles);


if strcmp(type, 'initialize')
    userData.statusM(userData.id).Visited = true;
    set(handles.figure1, 'UserData', userData)
end
