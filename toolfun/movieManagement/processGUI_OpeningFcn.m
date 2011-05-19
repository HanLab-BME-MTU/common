function processGUI_OpeningFcn(hObject, eventdata, handles, varargin)

[copyright openHelpFile] = userfcn_softwareConfig(handles);
set(handles.text_copyright, 'String', copyright)

userData = get(handles.figure1, 'UserData');
% Choose default command line output for thresholdProcessGUI
handles.output = hObject;

% Get main figure handle and process id
t = find(strcmp(varargin,'mainFig'));
userData.mainFig = varargin{t+1};
userData.procID = varargin{t+2};
userData.handles_main = guidata(userData.mainFig);

% Get current package and process
userData_main = get(userData.mainFig, 'UserData');
userData.crtPackage = userData_main.crtPackage;
userData.crtProc = userData.crtPackage.processes_{userData.procID};

% Get current process constructer
userData.procConstr = ...
    str2func(userData.crtPackage.processClassNames_{userData.procID});

% If process does not exist, create a default one in user data.
if isempty(userData.crtProc)
    userData.crtProc = userData.procConstr(userData_main.MD(userData_main.id), ...
        userData.crtPackage.outputDirectory_);
end

% Get icon infomation
userData.questIconData = userData_main.questIconData;
userData.colormap = userData_main.colormap;

% ---------------------- Channel Setup -------------------------

funParams = userData.crtProc.funParams_;

% Set up available input channels
set(handles.listbox_1, 'String', {userData_main.MD(userData_main.id).channels_.channelPath_},...
    'Userdata', 1: length(userData_main.MD(userData_main.id).channels_));

% Set up selected input data channels and channel index
parentI = find( userData.crtPackage.depMatrix_(userData.procID,:) );

if isempty(parentI) || ~isempty( userData.crtPackage.processes_{userData.procID} )
    
    % If process has no dependency, or process already exists, display saved channels
    set(handles.listbox_2, 'String', ...
        {userData_main.MD(userData_main.id).channels_(funParams.ChannelIndex).channelPath_}, ...
        'Userdata',funParams.ChannelIndex);
    
elseif isempty( userData.crtPackage.processes_{userData.procID} )
    % If new process
    empty = false;
    for i = parentI
        if isempty(userData.crtPackage.processes_{i})
            empty = true;
            break;
        end
    end
    
    if ~empty
        
        % If all dependent processes exist
        channelIndex = userData.crtPackage.processes_{parentI(1)}.funParams_.ChannelIndex;
        for i = 2: length(parentI)
            channelIndex = intersect(channelIndex, ...
                userData.crtPackage.processes_{parentI(i)}.funParams_.ChannelIndex);
        end
        
        if ~isempty(channelIndex)
            set(handles.listbox_2, 'String', ...
                {userData_main.MD(userData_main.id).channels(channelIndex).channelPath_}, ...
                'Userdata',channelIndex);
        end
    end
end
set(handles.checkbox_applytoall, 'Value', userData_main.applytoall(userData.procID));

% ----------------------Set up help icon------------------------

% Set up help icon
set(hObject,'colormap',userData.colormap);
% Set up package help. Package icon is tagged as '0'
axes(handles.axes_help);
Img = image(userData.questIconData);
set(gca, 'XLim',get(Img,'XData'),'YLim',get(Img,'YData'),...
    'visible','off','YDir','reverse');
set(Img,'ButtonDownFcn',@icon_ButtonDownFcn);
if openHelpFile
    set(Img, 'UserData', struct('class',class(userData.crtProc)))
end

% ----------------------------------------------------------------

% Update user data and GUI data
set(userData.mainFig, 'UserData', userData_main);
set(hObject, 'UserData', userData);

uicontrol(handles.pushbutton_done);
guidata(hObject, handles);

end