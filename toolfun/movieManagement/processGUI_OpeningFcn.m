function processGUI_OpeningFcn(hObject, eventdata, handles, string,varargin)
% Common initialization of concrete process GUIs
%
%       userData.mainFig - handle of main figure
%       userData.handles_main - 'handles' of main figure
%       userData.procID - The ID of process in the current package
%       userData.crtProc - handle of current process
%       userData.crtPackage - handles of current package
%       userData.procConstr - constructor of current process
%
%       userData.questIconData - help icon image information
%       userData.colormap - color map information
%
% Sebastien Besson May 2011

% Check input
ip = inputParser;
ip.addRequired('hObject',@ishandle);
ip.addRequired('eventdata',@(x) isstruct(x) || isempty(x));
ip.addRequired('handles',@isstruct);
ip.addRequired('string',@(x) isequal(x,'mainFig'));
ip.addOptional('mainFig',[],@ishandle);
ip.addOptional('procID',[],@isscalar);
ip.addOptional('process',[],@(x) isa(x,'function_handle'));
ip.parse(hObject,eventdata,handles,string,varargin{:});

% Retrieve userData and read function input 
userData = get(handles.figure1, 'UserData');
userData.mainFig=ip.Results.mainFig;
userData.procID = ip.Results.procID;
userData.procConstr=ip.Results.process;

% Set up copyright statement
[copyright openHelpFile] = userfcn_softwareConfig(handles);
set(handles.text_copyright, 'String', copyright)

% Get current package, movie data and process
userData.handles_main = guidata(userData.mainFig);
userData_main = get(userData.mainFig, 'UserData');
userData.MD = userData_main.MD(userData_main.id);
userData.crtPackage = userData_main.crtPackage;
userData.crtProc = userData.crtPackage.processes_{userData.procID};

% Get current process constructor
crtProcClassName = func2str(userData.procConstr);
crtProcName = eval([crtProcClassName '.getName']);
procString = [' Step ' num2str(userData.procID) ': ' crtProcName];
set(handles.text_processName,'String',procString);
figString = [' Setting - ' crtProcName];
set(handles.figure1,'Name',figString);

% If process does not exist, create a default one in user data.
if isempty(userData.crtProc)
    userData.crtProc = userData.procConstr(userData.MD, ...
        userData.crtPackage.outputDirectory_);
end

% Get icon infomation
userData.questIconData = userData_main.questIconData;
userData.colormap = userData_main.colormap;

% Check for multiple movies else
if numel(userData_main.MD) ==1
    set(handles.checkbox_applytoall,'Value',0,'Visible','off');
else
    set(handles.checkbox_applytoall, 'Value',...
        userData_main.applytoall(userData.procID));
end
uicontrol(handles.pushbutton_done);

% ----------------------Set up help icon------------------------

% Set up help icon
set(hObject,'colormap',userData.colormap);
% Set up package help. Package icon is tagged as '0'
set(handles.figure1,'CurrentAxes',handles.axes_help);
Img = image(userData.questIconData);
set(gca, 'XLim',get(Img,'XData'),'YLim',get(Img,'YData'),...
    'visible','off','YDir','reverse');
set(Img,'ButtonDownFcn',@icon_ButtonDownFcn);
if openHelpFile
    set(Img, 'UserData', struct('class',class(userData.crtProc)))
end

% ----------------------------------------------------------------

% Update user data and GUI data
set(hObject, 'UserData', userData);
end