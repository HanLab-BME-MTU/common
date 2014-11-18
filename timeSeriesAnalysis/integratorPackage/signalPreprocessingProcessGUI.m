function varargout = signalPreprocessingProcessGUI(varargin)
% signalPreprocessingProcessGUI M-file for signalPreprocessingProcessGUI.fig
%      signalPreprocessingProcessGUI, by itself, creates a new signalPreprocessingProcessGUI or raises the existing
%      singleton*.
%
%      H = signalPreprocessingProcessGUI returns the handle to a new signalPreprocessingProcessGUI or the handle to
%      the existing singleton*.
%
%      signalPreprocessingProcessGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in signalPreprocessingProcessGUI.M with the given input arguments.
%
%      signalPreprocessingProcessGUI('Property','Value',...) creates a new signalPreprocessingProcessGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before signalPreprocessingProcessGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to signalPreprocessingProcessGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help signalPreprocessingProcessGUI

% Last Modified by GUIDE v2.5 05-Mar-2012 17:00:21

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @signalPreprocessingProcessGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @signalPreprocessingProcessGUI_OutputFcn, ...
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


% --- Executes just before signalPreprocessingProcessGUI is made visible.
function signalPreprocessingProcessGUI_OpeningFcn(hObject,eventdata,handles,varargin)

processGUI_OpeningFcn(hObject, eventdata, handles, varargin{:});

userData=get(handles.figure1,'UserData');
if isempty(userData), userData = struct(); end
funParams = userData.crtProc.funParams_;

% Set up available movies
set(handles.listbox_availableMovies,'String',userData.ML.movieDataFile_, ...
    'UserData',1:numel(userData.ML.getMovies));

movieIndex = funParams.MovieIndex;

if ~isempty(movieIndex)
    movieString = userData.ML.movieDataFile_(movieIndex);
else
    movieString = {};
end

set(handles.listbox_selectedMovies,'String',movieString,...
    'UserData',movieIndex,'Callback',@(h,event)update_preview(h,event,guidata(h)));

signal=userData.ML.getMovies{1}.getSampledOutput;
nSignal=numel(signal);
% Create templates for GUI generation
allTrends=SignalPreprocessingProcess.getTrends;

createSignalText= @(pos,i) uicontrol(handles.uipanel_signal,'Style','text',...
    'Position',[20 pos 200 20],'Tag',['text_signal_' num2str(i)],...
    'String',signal(i).name,'HorizontalAlignment','left');
createOutlierBox= @(pos,i) uicontrol(handles.uipanel_signal,'Style','edit',...
    'Position',[250 pos 60 20],'Tag',['edit_kSigma_' num2str(i)],...
    'BackgroundColor','w','String',funParams.kSigma(i));
createTrendMenu= @(pos,i) uicontrol(handles.uipanel_signal,'Style','popupmenu',...
    'Position',[330 pos 140 20],'Tag',['popupmenu_trendType_' num2str(i)],...
    'String',{allTrends.name},'UserData',[allTrends.type],...
    'Value',find([allTrends.type]==funParams.trendType(i)));


% Create checkboxes for samplable processes
hPosition =10;
for i = nSignal:-1:1
    createSignalText(hPosition,i);
    createTrendMenu(hPosition,i);
    createOutlierBox(hPosition,i);
    hPosition=hPosition+30;
end
userData.previewFig=-ones(nSignal,1);


% Add cosmetic title for signal
uicontrol(handles.uipanel_signal,'Style','text',...
    'Position',[20 hPosition 150 30],'String','Signal name');
uicontrol(handles.uipanel_signal,'Style','text',...
    'Position',[250 hPosition 60 30],'String','Outlier threshold');
uicontrol(handles.uipanel_signal,'Style','text',...
    'Position',[330 hPosition 140 30],'String','Trend to remove');

% Fix GUI position/size
a=get(get(handles.uipanel_signal,'Children'),'Position');
P=vertcat(a{:});
panelSize = [max(P(:,1)+P(:,3))+10 max(P(:,2)+P(:,4))+20];
pos = get(handles.uipanel_signal,'Position');
dh= panelSize(2) - pos(4);
set(handles.figure1,'Position',get(handles.figure1,'Position')+[0 -dh 0 dh]);
set(handles.uipanel_signal,'Position',get(handles.uipanel_signal,'Position')+[0 0 0 dh]);

uicontrols = {'uipanel_movies','axes_help','text_processName','text_copyright'};
for i=1:numel(uicontrols)
    set(handles.(uicontrols{i}),'Position',get(handles.(uicontrols{i}),'Position')+[0 dh 0 0]);
end

% Set windows selection parameters
set(handles.edit_BandMin,'String',funParams.BandMin);
set(handles.edit_BandMax,'String',funParams.BandMax);
userData.SliceIndex=funParams.SliceIndex;

% Set energy calculation parameters
set(handles.edit_nBoot,'String',funParams.nBoot);
set(handles.edit_alpha,'String',funParams.alpha);


% Update handles structure and attach it to the main figure
handles = guihandles(handles.figure1);

% Choose default command line output for signalPreprocessingProcessGUI
handles.output = hObject;

% Update user data and GUI data 
set(hObject, 'UserData', userData);
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = signalPreprocessingProcessGUI_OutputFcn(~, ~, handles) 
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
if isempty(userData), userData = struct(); end

delete(userData.helpFig(ishandle(userData.helpFig)));
delete(userData.previewFig(ishandle(userData.previewFig)));


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

% --- Executes on button press in checkbox_all.
function checkbox_all_Callback(hObject, eventdata, handles)

% Retrieve available channels properties
props = get(handles.listbox_availableMovies, {'String','UserData'});
if isempty(props{1}), return; end

% Update selected channels
if get(hObject,'Value')
    set(handles.listbox_selectedMovies, 'String', props{1},...
        'UserData',props{2});
else
    set(handles.listbox_selectedMovies, 'String', {}, 'UserData',[], 'Value',1);
end

% --- Executes on button press in pushbutton_selectMovies.
function pushbutton_select_Callback(hObject, eventdata, handles)


% Retrieve  Movies properties
availableProps = get(handles.listbox_availableMovies, {'String','UserData','Value'});
selectedProps = get(handles.listbox_selectedMovies, {'String','UserData'});

% Find new elements and set them to the selected listbox
newID = availableProps{3}(~ismember(availableProps{1}(availableProps{3}),selectedProps{1}));
selectedMovies = horzcat(selectedProps{1}',availableProps{1}(newID)');
selectedData = horzcat(selectedProps{2}, availableProps{2}(newID));
set(handles.listbox_selectedMovies, 'String', selectedMovies, 'UserData', selectedData);


% --- Executes on button press in pushbutton_deleteMovies.
function pushbutton_delete_Callback(hObject, eventdata, handles)

% Get selected properties and returin if empty
selectedProps = get(handles.listbox_selectedMovies, {'String','UserData','Value'});
if isempty(selectedProps{1}) || isempty(selectedProps{3}),return; end

% Delete selected item
selectedProps{1}(selectedProps{3}) = [ ];
selectedProps{2}(selectedProps{3}) = [ ];
set(handles.listbox_selectedMovies, 'String', selectedProps{1},'UserData',selectedProps{2},...
    'Value',max(1,min(selectedProps{3},numel(selectedProps{1}))));

% --- Executes on button press in pushbutton_done.
function pushbutton_done_Callback(hObject, eventdata, handles)


% Check user input
if isempty(get(handles.listbox_selectedMovies, 'String'))
    errordlg('Please select at least one input process from ''Available Movies''.','Setting Error','modal')
    return;
end
funParams.MovieIndex = get(handles.listbox_selectedMovies, 'UserData');

h=findobj('-regexp','Tag','text_signal_(\d)');
nSignal=numel(h);
kSigma=NaN(nSignal,1);
for i=1:nSignal
    kSigma(i) = str2double(get(handles.(['edit_kSigma_'  num2str(i)]),'String'));
    if isnan(kSigma(i)) || ~ismember(kSigma(i),1:10)
        errordlg(['Please enter a valid outlie threshold for the '...
            get(handles.(['text_signal_' num2str(i)]),'String')],...
            'Setting error','modal');
        return;
    end
end
funParams.kSigma = kSigma;

% Retrieve trendType
funParams.trendType=NaN(nSignal,1);
for i=1:nSignal
    props = get(handles.(['popupmenu_trendType_' num2str(i)]),{'UserData','Value'});
    funParams.trendType(i)=props{1}(props{2});
end

nBoot = str2double(get(handles.edit_nBoot,'String'));
if isnan(nBoot) || nBoot<=0 || round(nBoot) ~= nBoot
    errordlg(['Please enter a valid value for the ' get(handles.text_nBoot,'String')],...
        'Setting error','modal');
    return;
end
funParams.nBoot = nBoot;

alpha = str2double(get(handles.edit_alpha,'String'));
if isnan(alpha) || alpha<=0 || alpha >=1
    errordlg(['Please enter a valid value for the ' get(handles.text_alpha,'String')],...
        'Setting error','modal');
    return;
end
funParams.alpha = alpha;




bandMin = str2double(get(handles.edit_BandMin,'String'));
if isnan(bandMin) || bandMin<1
    errordlg('Please enter a valid value for the minimum band of windows to correlate',...
        'Setting error','modal');
    return;
end

bandMax = str2double(get(handles.edit_BandMax,'String'));
if isnan(bandMax) 
    errordlg('Please enter a valid value for the maximum band of windows to correlate',...
        'Setting error','modal');
    return;
end

funParams.MovieIndex = get(handles.listbox_selectedMovies, 'UserData');
funParams.BandMin=bandMin;
funParams.BandMax=bandMax;

userData = get(handles.figure1, 'UserData');
if ishandle(userData.previewFig)
    alphamask = get(findobj(userData.previewFig(1),'Tag','sliceIndexImage'),'AlphaData');
    userData.SliceIndex{userData.movieID}=logical(sum(alphamask,2));
    delete(userData.previewFig(1))
end
funParams.SliceIndex=userData.SliceIndex;

% Set parameters
processGUI_ApplyFcn(hObject, eventdata, handles,funParams);


% --- Executes on button press in checkbox_selectSlices.
function update_preview(hObject, eventdata, handles)

userData=get(handles.figure1,'UserData');
if isempty(userData), userData = struct(); end

% Delete figure if checkbox is unselected
if ~get(handles.checkbox_selectSlices,'Value')
    delete(userData.previewFig(ishandle(userData.previewFig)));
    return;
end

if isequal(hObject,handles.listbox_selectedMovies) && ishandle(userData.previewFig(1)), 
    alphamask = get(findobj(userData.previewFig(1),'Tag','sliceIndexImage'),'AlphaData');
    userData.SliceIndex{userData.movieID}=logical(sum(alphamask,2));
end

% Retrieve new movie ID
movieProps = get(handles.listbox_selectedMovies,{'UserData','Value'});
userData.movieID=movieProps{1}(movieProps{2});
movie = userData.ML.getMovies{userData.movieID};

% Retrieve input 
signal = movie.getSampledOutput;
nSignal = numel(signal);

% Create a mask using the slice index
alphamask = repmat(userData.SliceIndex{userData.movieID},1,movie.nFrames_);
hAxes = -ones(nSignal,1);
h = -ones(nSignal,1);
userData_fig.mainFig=handles.figure1;
userData_fig.alphamask = alphamask;

for i=find(~ishandle(userData.previewFig))'
    userData.previewFig(i)=figure('NumberTitle','off','Name',signal(i).name);
end

for i=1:nSignal
    figure(userData.previewFig(i));
    clf(userData.previewFig(i));
    hAxes(i) = axes('HitTest','on','ButtonDownFcn',@(h,event)editSliceIndex(h,event,guidata(h)));
    procID=signal(i).processIndex;
    if ~isempty(signal(i).channelIndex), 
        drawArgs={signal(i).channelIndex,signal(i).outputIndex}; 
    else
        drawArgs={};
    end
    h(i) = movie.processes_{procID}.draw(drawArgs{:});  
    set(h(i),'Tag','sliceIndexImage',...
        'HitTest','off','AlphaData',alphamask,'AlphaDataMapping','none');
    set(userData.previewFig(i),'UserData',userData_fig,'DeleteFcn',@(h,event)closeGraphFigure(h,event));
end
linkaxes(hAxes);

set(handles.figure1, 'UserData', userData);

function editSliceIndex(src,eventdata,handles)

% Retrieve window index corresponding to point poistion
point = get(src,'CurrentPoint');
windowIndex=floor(point(1,2));

f =  get(src,'Parent');
userData = get(f,'UserData');
userData.windowStart=windowIndex;
userData.alphavalue = ~userData.alphamask(windowIndex,1);

% userData.alphamask(windowIndex,:)=.2;
h=findobj(f,'Tag','sliceIndexImage');
set(h,'AlphaData',userData.alphamask);
set(f,'UserData',userData,'WindowButtonMotionFcn',@updateAlphaMask,...
    'WindowButtonUpFcn',@saveAlphaMask);
   
function saveAlphaMask(src,event)

userData = updateAlphaMask(src);
if isempty(userData), userData = struct(); end
% Save data in the figure
set(src,'UserData',userData);
set(src,'WindowButtonMotionFcn',[],'WindowButtonUpFcn',[]);
 
function userData=updateAlphaMask(src,event)

% Retrieve current point
userData = get(src,'UserData');
if isempty(userData), userData = struct(); end
point = get(gca,'CurrentPoint');
windowIndex=floor(point(1,2));

% Scale within the axes limits and the mask range
yLim=floor(get(gca,'YLim'));
windowIndex=max(min(yLim(2),windowIndex),yLim(1));
windowIndex=max(min(size(userData.alphamask,1),windowIndex),1);

% Update selected range
if windowIndex>=userData.windowStart
    windowRange=userData.windowStart:windowIndex;
else
    windowRange=windowIndex:userData.windowStart;
end

% Update graphic alpha mask
userData = get(src,'UserData');
userData.alphamask(windowRange,:)=userData.alphavalue;
set(findobj(0,'Tag','sliceIndexImage'),'AlphaData',userData.alphamask);

function closeGraphFigure(src,event)

% Update SliceIndex parameters in main process figure
userData = get(src,'UserData');
if isempty(userData), userData = struct(); end
userData_main =get(userData.mainFig,'UserData');
handles_main = guidata(userData.mainFig);
userData_main.SliceIndex{userData_main.movieID}=logical(sum(userData.alphamask,2));

% Update main figure
delete(userData_main.previewFig(ishandle(userData_main.previewFig)));
set(userData.mainFig,'UserData',userData_main);

set(handles_main.checkbox_selectSlices,'Value',0);   
