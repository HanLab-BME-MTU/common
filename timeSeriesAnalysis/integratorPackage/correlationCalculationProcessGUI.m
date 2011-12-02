function varargout = correlationCalculationProcessGUI(varargin)
% correlationCalculationProcessGUI M-file for correlationCalculationProcessGUI.fig
%      correlationCalculationProcessGUI, by itself, creates a new correlationCalculationProcessGUI or raises the existing
%      singleton*.
%
%      H = correlationCalculationProcessGUI returns the handle to a new correlationCalculationProcessGUI or the handle to
%      the existing singleton*.
%
%      correlationCalculationProcessGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in correlationCalculationProcessGUI.M with the given input arguments.
%
%      correlationCalculationProcessGUI('Property','Value',...) creates a new correlationCalculationProcessGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before correlationCalculationProcessGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to correlationCalculationProcessGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help correlationCalculationProcessGUI

% Last Modified by GUIDE v2.5 18-Nov-2011 17:06:46

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @correlationCalculationProcessGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @correlationCalculationProcessGUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before correlationCalculationProcessGUI is made visible.
function correlationCalculationProcessGUI_OpeningFcn(hObject,eventdata,handles,varargin)

processGUI_OpeningFcn(hObject, eventdata, handles, varargin{:});

userData=get(handles.figure1,'UserData');
funParams = userData.crtProc.funParams_;

% Set up available input channels
set(handles.listbox_availableMovies,'String',userData.ML.movieDataFile_, ...
    'UserData',1:numel(userData.ML.movies_));

movieIndex = funParams.MovieIndex;

if ~isempty(movieIndex)
    movieString = userData.ML.movieDataFile_(movieIndex);
else
    movieString = {};
end

set(handles.listbox_selectedMovies,'String',movieString,...
    'UserData',movieIndex,'Callback',@(h,event)update_preview(h,event,guidata(h)));

% Set up available input processes
allProc = userData.crtProc.getTimeSeriesProcesses();
allProcString = cellfun(@(x) eval([x '.getName']),allProc,'UniformOutput',false);
set(handles.listbox_availableProcesses,'String',allProcString,'UserData',allProc);

% Set up selected input processes
selProc = funParams.ProcessName;
selProcString = cellfun(@(x) eval([x '.getName']),selProc,'UniformOutput',false);
set(handles.listbox_selectedProcesses,'String',selProcString,'UserData',selProc);

set(handles.edit_BandMin,'String',funParams.BandMin);
set(handles.edit_BandMax,'String',funParams.BandMax);
userData.SliceIndex=funParams.SliceIndex;
userData.previewFig=-1;

% Choose default command line output for correlationCalculationProcessGUI
handles.output = hObject;

% Update user data and GUI data
set(hObject, 'UserData', userData);
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = correlationCalculationProcessGUI_OutputFcn(~, ~, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in pushbutton_cancel.
function pushbutton_cancel_Callback(~, ~, handles)
% Delete figure
delete(handles.figure1);

% --- Executes during object deletion, before destroying properties.
function figure1_DeleteFcn(hObject, ~, handles)
% Notify the package GUI that the setting panel is closed
userData = get(handles.figure1, 'UserData');

if ishandle(userData.helpFig), delete(userData.helpFig); end
if ishandle(userData.previewFig), delete(userData.previewFig); end

set(handles.figure1, 'UserData', userData);
guidata(hObject,handles);


% --- Executes on key press with focus on pushbutton_done and none of its controls.
function pushbutton_done_KeyPressFcn(~, eventdata, handles)

if strcmp(eventdata.Key, 'return')
    pushbutton_done_Callback(handles.pushbutton_done, [], handles);
end

% --- Executes on key press with focus on figure1 and none of its controls.
function figure1_KeyPressFcn(~, eventdata, handles)

if strcmp(eventdata.Key, 'return')
    pushbutton_done_Callback(handles.pushbutton_done, [], handles);
end

% --- Executes on button press in checkbox_allMovies.
function checkbox_all_Callback(hObject, eventdata, handles)

% Identify listbox and retrieve handles
tokens = regexp(get(hObject,'Tag'),'^checkbox_all(.*)$','tokens');
listbox_available= handles.(['listbox_available' tokens{1}{1}]);
listbox_selected= handles.(['listbox_selected' tokens{1}{1}]);

% Retrieve available properties
availableProps = get(listbox_available, {'String','UserData'});
if isempty(availableProps{1}), return; end

if get(hObject,'Value')
    set(listbox_selected, 'String', availableProps{1},'UserData',availableProps{2});
else
    set(listbox_selected, 'String', {}, 'UserData',[], 'Value',1);
end

% --- Executes on button press in pushbutton_selectMovies.
function pushbutton_select_Callback(hObject, eventdata, handles)

% Identify listbox and retrieve handles
tokens = regexp(get(hObject,'Tag'),'^pushbutton_select(.*)$','tokens');
listbox_available= handles.(['listbox_available' tokens{1}{1}]);
listbox_selected= handles.(['listbox_selected' tokens{1}{1}]);

% Get handles properties
availableProps = get(listbox_available, {'String','UserData'});
selectedProps = get(listbox_selected, {'String','UserData'});
ID = get(listbox_available, 'Value');

% Update selected listbox properties
newChanID = ID(~ismember(availableProps{1}(ID),selectedProps{1}));
selectedString = vertcat(selectedProps{1},availableProps{1}(newChanID));
selectedData = horzcat(selectedProps{2}, availableProps{2}(newChanID));

set(listbox_selected, 'String', selectedString, 'Userdata', selectedData);


% --- Executes on button press in pushbutton_deleteMovies.
function pushbutton_delete_Callback(hObject, eventdata, handles)

% Identify listbox and retrieve handles
tokens = regexp(get(hObject,'Tag'),'^pushbutton_delete(.*)$','tokens');
listbox_selected= handles.(['listbox_selected' tokens{1}{1}]);

% Get selected properties and returin if empty
selectedProps = get(listbox_selected, {'String','UserData','Value'});
if isempty(selectedProps{1}) || isempty(selectedProps{3}),return; end

% Delete selected item
selectedProps{1}(selectedProps{3}) = [ ];
selectedProps{2}(selectedProps{3}) = [ ];
set(listbox_selected, 'String', selectedProps{1},'UserData',selectedProps{2},...
    'Value',max(1,min(selectedProps{3},numel(selectedProps{1}))));
update_preview(hObject, eventdata, handles)

% --- Executes on button press in checkbox_selectSlices.
function update_preview(hObject, eventdata, handles)

userData=get(handles.figure1,'UserData');

% Delete figure if checkbox is unselected
if ~get(handles.checkbox_selectSlices,'Value') 
    if ishandle(userData.previewFig), delete(userData.previewFig); end
    return; 
end

if isequal(hObject,handles.checkbox_selectSlices)
    % Create new figure
    userData.previewFig=figure('NumberTitle','off','Name','Select/unselect the windows of interest');
else
    % Save slice index of current movie and clear figure
    figure(userData.previewFig);
    userData_fig=get(userData.previewFig,'UserData');
    userData.SliceIndex{userData.movieID}=logical(sum(userData_fig.alphamask,2));
    clf(userData.previewFig)
end


% Retrieve movie id and process names
movieProps = get(handles.listbox_selectedMovies,{'UserData','Value'});
userData.movieID=movieProps{1}(movieProps{2});
p.ProcessName = get(handles.listbox_selectedProcesses,'UserData');

% Retrieve input 
corrProc = CorrelationCalculationProcess(userData.ML.movies_{userData.movieID},'');
parseProcessParams(corrProc,p);
input = corrProc.getInput;
nInput = numel(input);


% Create a mask using the slice index
alphamask = repmat(userData.SliceIndex{userData.movieID},1,...
    userData.ML.movies_{userData.movieID}.nFrames_);
hAxes = -ones(nInput,1);
h = -ones(nInput,1);
userData_fig.mainFig=handles.figure1;
userData_fig.alphamask = alphamask;
userData_fig.crtAxes=[];
for i=1:nInput;
    hAxes(i) = axes('Position',[(i-1)/nInput 0.1 1/nInput .8],'HitTest','on',...
        'ButtonDownFcn',@(h,event)editSliceIndex(h,event,guidata(h)));
    axesArgs={'hAxes',hAxes(i)};
    procID=input(i).processIndex;
    if ~isempty(input(i).channelIndex)
        chanArgs={input(i).channelIndex};
    else
        chanArgs={};
    end
    h(i) = userData.ML.movies_{userData.movieID}.processes_{procID}.draw(chanArgs{:},axesArgs{:});
    hCbar = findobj(userData.previewFig,'Tag','Colorbar');
    delete(hCbar);    
    set(h(i),'HitTest','off','AlphaData',alphamask,'AlphaDataMapping','none');
end
linkaxes(hAxes);

set(userData.previewFig,'UserData',userData_fig,....
    'DeleteFcn',@(h,event)closeGraphFigure(h,event));
set(handles.figure1, 'UserData', userData);


function editSliceIndex(src,eventdata,handles)


point = get(src,'CurrentPoint');
windowIndex=ceil(point(1,2));

f =  get(src,'Parent');
userData = get(f,'UserData');
userData.crtAxes=src;
userData.windowStart=windowIndex;
userData.alphavalue = ~userData.alphamask(windowIndex,1);

% userData.alphamask(windowIndex,:)=.2;
h=findobj(f,'Type','Image');
set(h,'AlphaData',userData.alphamask);
set(f,'UserData',userData,'WindowButtonMotionFcn',@moveSliceIndex,...
    'WindowButtonUpFcn',@releaseSliceIndex);
   

function moveSliceIndex(src,event)

% Retrieve current point
userData = get(src,'UserData');
point = get(userData.crtAxes,'CurrentPoint');
windowIndex=ceil(point(1,2));

% Scale within the axes limits and the mask range
yLim=floor(get(userData.crtAxes','YLim'));
windowIndex=max(min(yLim(2),windowIndex),yLim(1));
windowIndex=max(min(size(userData.alphamask,1),windowIndex),1);

% Update selected range
if windowIndex>=userData.windowStart
    windowRange=userData.windowStart:windowIndex;
else
    windowRange=windowIndex:userData.windowStart;
end
 
% Update graphic alpha mask
userData.alphamask(windowRange,:)=userData.alphavalue;
set(findobj(src,'Type','Image'),'AlphaData',userData.alphamask);

function releaseSliceIndex(src,event)

% Retrieve current point
userData = get(src,'UserData');
point = get(userData.crtAxes,'CurrentPoint');
windowIndex=floor(point(1,2));

% Scale within the axes limits and the mask range
yLim=floor(get(userData.crtAxes','YLim'));
windowIndex=max(min(yLim(2),windowIndex),yLim(1));
windowIndex=max(min(size(userData.alphamask,1),windowIndex),1);

% Update selected range
if windowIndex>=userData.windowStart
    windowRange=userData.windowStart:windowIndex;
else
    windowRange=windowIndex:userData.windowStart;
end
 
% Update graphic alpha mask
userData.alphamask(windowRange,:)=userData.alphavalue;
set(findobj(src,'Type','Image'),'AlphaData',userData.alphamask);

% Save data in the figure
set(src,'UserData',userData);
set(src,'WindowButtonMotionFcn',[],'WindowButtonUpFcn',[]);
 
function closeGraphFigure(src,event)

% Update SliceIndex parameters in main process figure
userData = get(src,'UserData');
userData_main =get(userData.mainFig,'UserData');
handles_main = guidata(userData.mainFig);
userData_main.SliceIndex{userData_main.movieID}=logical(sum(userData.alphamask,2));
set(userData.mainFig,'UserData',userData_main);
set(handles_main.checkbox_selectSlices,'Value',0);   


% --- Executes on button press in pushbutton_done.
function pushbutton_done_Callback(hObject, eventdata, handles)

% Check user input
if isempty(get(handles.listbox_selectedMovies, 'String'))
    errordlg('Please select at least one input process from ''Available Movies''.','Setting Error','modal')
    return;
end

if isempty(get(handles.listbox_selectedProcesses, 'String'))
    errordlg('Please select at least one input process from ''Available Processes''.','Setting Error','modal')
    return;
end

bandMin = str2double(get(handles.edit_BandMin,'String'));
if isnan(bandMin) || bandMin<1
    errordlg('Please enter a valid value for the layer of cells to be used',...
        'Setting error','modal');
    return;
end

bandMax = str2double(get(handles.edit_BandMax,'String'));
if isnan(bandMax) 
    errordlg('Please enter a valid value for the maximum layer of cells to be used',...
        'Setting error','modal');
    return;
end

funParams.MovieIndex = get(handles.listbox_selectedMovies, 'UserData');
funParams.ProcessName = get(handles.listbox_selectedProcesses, 'UserData');
funParams.BandMin=bandMin;
funParams.BandMax=bandMax;

userData = get(handles.figure1, 'UserData');
if ishandle(userData.previewFig)
    userData_fig=get(userData.previewFig,'UserData');
    userData.SliceIndex{userData.movieID}=logical(sum(userData_fig.alphamask,2));
    delete(userData.previewFig)
end
funParams.SliceIndex=userData.SliceIndex;

% Process Sanity check ( only check underlying data )
userData = get(handles.figure1, 'UserData');
try
    userData.crtProc.sanityCheck;
catch ME

    errordlg([ME.message 'Please double check your data.'],...
                'Setting Error','modal');
    return;
end

% Set parameters
processGUI_ApplyFcn(hObject, eventdata, handles,funParams);
