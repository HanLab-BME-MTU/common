function varargout = trackingVisualGUI(varargin)
% TRACKINGVISUALGUI M-file for trackingVisualGUI.fig
%      TRACKINGVISUALGUI, by itself, creates a new TRACKINGVISUALGUI or raises the existing
%      singleton*.
%
%      H = TRACKINGVISUALGUI returns the handle to a new TRACKINGVISUALGUI or the handle to
%      the existing singleton*.
%
%      TRACKINGVISUALGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TRACKINGVISUALGUI.M with the given input arguments.
%
%      TRACKINGVISUALGUI('Property','Value',...) creates a new TRACKINGVISUALGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before trackingVisualGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to trackingVisualGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help trackingVisualGUI

% Last Modified by GUIDE v2.5 23-Mar-2011 11:15:41

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @trackingVisualGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @trackingVisualGUI_OutputFcn, ...
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


% --- Executes just before trackingVisualGUI is made visible.
function trackingVisualGUI_OpeningFcn(hObject, eventdata, handles, varargin)

% userData.showDetFig = trackingVisualGUI('mainFig',handles.figure1, procID);
%
% Available tools 
% UserData data:
%       userData.MD - 1x1 the current movie data
%       userData.mainFig - handle of main figure
%       userData.handles_main - 'handles' of main figure
%       userData.procID - The ID of process in the current package
%       userData.crtProc - handle of current process
%       
%       userData.chan - the channel index to display
%       userData.firstframe - the first frame who has tracking result
%       userData.lastframe - the last frame who has tracking result
%       userData.file - maximum frame number
%

[copyright openHelpFile] = userfcn_softwareConfig(handles);
set(handles.text_copyright, 'String', copyright)

userData = get(handles.figure1, 'UserData');
% Choose default command line output for detectionVisualGUI
handles.output = hObject;

% Get main figure handle and process id
t = find(strcmp(varargin,'mainFig'));
userData.mainFig = varargin{t+1};
userData.procID = varargin{t+2};
userData.handles_main = guidata(userData.mainFig);
userData.userData_main = get(userData.handles_main.figure1, 'UserData');

% Get current package and process
userData_main = get(userData.mainFig, 'UserData');
userData.MD = userData_main.MD(userData_main.id);  % Get the current Movie Data
userData.crtPackage = userData_main.crtPackage;
userData.crtProc = userData.crtPackage.processes_{userData.procID};

% Get channel index
chan = [];
for i = 1:length(userData.MD.channels_)
    if userData.crtProc.checkChannelOutput(i)
        chan = i; 
        break
    end
end

assert(~isempty(chan), 'User-defined: the process does not have output yet.')
userData.chan = chan;

% Make sure detection output is valid
load(userData.crtProc.outFilePaths_{chan},'tracksFinal');
if isempty(tracksFinal)
   error('User-defined: there is no detection information in the output variable.') 
end

userData.tracksFinal = tracksFinal;

allEvents = vertcat(userData.tracksFinal.seqOfEvents);
userData.firstframe = min(allEvents(:,1));
userData.lastframe = max(allEvents(:,1));

userData.toolName = {'plotTracks2D', 'plotCompTrack', 'overlayTracksMovieNew'};

% ------------------ Update input parameters ------------------------

visualParams = userData.crtProc.visualParams_;

% Tool 1: plotTrakcs2D
visualParams.pt2D.timeRange = [userData.firstframe, userData.lastframe]; 

% Tool 2: plotCompTrack

% Tool 3: overlayTracksMovieNew
visualParams.otmn.startend = [userData.firstframe, userData.lastframe]; 

file = userData.MD.getImageFileNames(chan);
visualParams.otmn.firstImageFile = [userData.MD.channels_(chan).channelPath_ filesep file{1}{1}];

if any(strcmp(visualParams.otmn.dir2saveMovie, {userData.MD.channels_(:).channelPath_}))
    visualParams.otmn.dir2saveMovie = userData.MD.channels_(chan).channelPath_;
end


userData.crtProc.setVisualParams(visualParams)

% ------------------- GUI Notes set-up -------------------------------

% text_movie
str = userData.MD.channels_(chan).channelPath_;
limit = 180;
if length(str)>limit
    str = ['... ' str(end-limit:end)];
end
set(handles.text_movie, 'String', str)

% text_framenum
set(handles.text1_framenum, 'String', ['( Tracks availabe from frame ',num2str(userData.firstframe),' to ',num2str(userData.lastframe),' )'])
set(handles.text3_framenum, 'String', ['( Tracks availabe from frame ',num2str(userData.firstframe),' to ',num2str(userData.lastframe),' )'])
set(handles.text2_tracknum, 'String', ['out of ',num2str(length(userData.tracksFinal)),' tracks'])

% popupmenu set-up
colorStr = {'Color-code Time (G->B->R)', 'Rotate Through 7 Colors', 'Rotate Through 23 Colors','Black', 'Blue', 'Red'};
colorUserData = {'1','2','3','k','b','r'};
set(handles.popupmenu1_color, 'String', colorStr, 'UserData', colorUserData)

markerStr = {'Plus (+)', 'Circle (o)', 'Asterisk (*)','Point (.)', 'Cross (x)', 'Square', 'Diamond', 'none'};
markerUserData = {'+','o','*','.','x','square', 'diamond', 'none'};
set(handles.popupmenu1_marker, 'String', markerStr, 'UserData', markerUserData)



% -------------------- Set Parameters -------------------------------

visualParams = userData.crtProc.visualParams_;

% Tool 1: plotTrakcs2D

set(handles.edit1_min, 'String', num2str(visualParams.pt2D.timeRange(1)))
set(handles.edit1_max, 'String', num2str(visualParams.pt2D.timeRange(2)))

set(handles.popupmenu1_color, 'Value', find(strcmp(colorUserData, visualParams.pt2D.colorTime)))
set(handles.popupmenu1_marker, 'Value', find(strcmp(markerUserData, visualParams.pt2D.markerType)))

set(handles.checkbox1_indicateSE, 'Value', visualParams.pt2D.indicateSE)
set(handles.checkbox1_flipXY, 'Value', visualParams.pt2D.flipXY)
set(handles.checkbox1_ask4sel, 'Value', visualParams.pt2D.ask4sel)

set(handles.edit1_offset1, 'String', num2str(visualParams.pt2D.offset(1)))
set(handles.edit1_offset2, 'String', num2str(visualParams.pt2D.offset(2)))

if visualParams.pt2D.newFigure
    set(handles.edit1_image, 'String', visualParams.pt2D.imageDir)
end


% Tool 2: plotCompTrack

set(handles.edit2_trackid, 'String', num2str(visualParams.pct.trackid))
set(handles.checkbox2_plotX, 'Value', visualParams.pct.plotX)
set(handles.checkbox2_plotY, 'Value', visualParams.pct.plotY)
set(handles.checkbox2_plotA, 'Value', visualParams.pct.plotA)
set(handles.checkbox2_inOneFigure, 'Value', visualParams.pct.inOneFigure)


% Tool 3: overlayTracksMovieNew

set(handles.edit3_min, 'String', num2str(visualParams.otmn.startend(1)))
set(handles.edit3_max, 'String', num2str(visualParams.otmn.startend(2)))
set(handles.edit3_dragtailLength, 'String', num2str(visualParams.otmn.dragtailLength))
set(handles.edit3_filterSigma, 'String', num2str(visualParams.otmn.filterSigma))

% Set default movie name as the name of result MAT file
% if isempty(visualParams.otmn.movieName)

     fName = userData.crtProc.funParams_.saveResults.filename;
     [x1 x2 x3 x4] = getFilenameBody(fName);
     visualParams.otmn.movieName = [x2 x3];
     userData.crtProc.setVisualParams(visualParams)
     visualParams = userData.crtProc.visualParams_;
% end

if ~visualParams.otmn.saveMovie
    
    set(handles.checkbox3_save, 'Value', 0)
    set(handles.text3_filename, 'Enable', 'off')
    set(handles.edit3_filename, 'Enable', 'off', 'String', '')
    set(handles.text3_mov, 'Enable', 'off')
    set(handles.pushbutton3_path, 'Enable', 'off')
    set(handles.edit3_path, 'Enable', 'off')
    
else
    set(handles.checkbox3_save, 'Value', 1)
    set(handles.edit3_filename, 'String', visualParams.otmn.movieName)
    set(handles.edit3_path, 'String', visualParams.otmn.dir2saveMovie)
    
end

set(handles.checkbox3_onlyTracks, 'Value',  visualParams.otmn.onlyTracks)
set(handles.checkbox3_colorTracks, 'Value',  visualParams.otmn.colorTracks)
set(handles.checkbox3_classifyGaps, 'Value',  visualParams.otmn.classifyGaps)
set(handles.checkbox3_highlightES, 'Value',  visualParams.otmn.highlightES)
set(handles.checkbox3_classifyLft, 'Value',  visualParams.otmn.classifyLft)

if visualParams.otmn.onlyTracks || visualParams.otmn.colorTracks
    
    set(handles.checkbox3_classifyGaps, 'Value', 0, 'Enable', 'off')
    set(handles.checkbox3_highlightES, 'Value', 0, 'Enable', 'off')
    set(handles.checkbox3_classifyLft, 'Value', 0, 'Enable', 'off')
    set(handles.text44, 'Enable', 'off')    
end

if ~isempty(visualParams.otmn.imageRange)
    
    set(handles.edit3_minx, 'String', num2str(visualParams.otmn.imageRange(1,1)))
    set(handles.edit3_maxx, 'String', num2str(visualParams.otmn.imageRange(1,2)))
    set(handles.edit3_miny, 'String', num2str(visualParams.otmn.imageRange(2,1)))
    set(handles.edit3_maxy, 'String', num2str(visualParams.otmn.imageRange(2,2)))
end

switch visualParams.otmn.showRaw
    
    case 0 
        set(handles.radiobutton3_display_0, 'Value', 1)
    case 1
        set(handles.radiobutton3_display_1, 'Value', 1)
    case 2
        set(handles.radiobutton3_display_2, 'Value', 1)
    otherwise
        error('User-defined: the parameter is incorrect.')
end

switch visualParams.otmn.intensityScale
    
    case 0 
        set(handles.radiobutton3_scale_0, 'Value', 1)
    case 1
        set(handles.radiobutton3_scale_1, 'Value', 1)
    case 2
        set(handles.radiobutton3_scale_2, 'Value', 1)
    otherwise
        error('User-defined: the parameter is incorrect.')
end

% Get icon infomation
userData.questIconData = userData.userData_main.questIconData;
userData.colormap = userData.userData_main.colormap;

% ----------------------Set up help icon------------------------    

% Set up help icon
set(hObject,'colormap',userData.colormap);

for i = 1:length(userData.toolName)

    % Set up package help. Package icon is tagged as '0'
    eval (['axes(handles.axes_help_' num2str(i) ')'])
    Img = image(userData.questIconData); 
    set(gca, 'XLim',get(Img,'XData'),'YLim',get(Img,'YData'),...
        'visible','off','YDir','reverse');
    set(Img,'ButtonDownFcn',@icon_ButtonDownFcn);
    
    if openHelpFile
        set(Img, 'UserData', struct('class', userData.toolName{i}))
    else
        set(Img, 'UserData', 'Please refer to help file.')
    end

end

% Update user data and GUI data
set(hObject, 'UserData', userData);
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = trackingVisualGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton_done.
function pushbutton_done_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_done (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

delete(handles.figure1)


% --- Executes on button press in pushbutton3_display.
function pushbutton3_display_Callback(hObject, eventdata, handles)

userData = get(handles.figure1, 'UserData');

min = get(handles.edit3_min, 'String');
max = get(handles.edit3_max, 'String');
filterSigma = get(handles.edit3_filterSigma, 'String');
dragtailLength = get(handles.edit3_dragtailLength, 'String');
filename = get(handles.edit3_filename, 'String');
path = get(handles.edit3_path, 'String');

minx = get(handles.edit3_minx, 'String');
maxx = get(handles.edit3_maxx, 'String');
miny = get(handles.edit3_miny, 'String');
maxy = get(handles.edit3_maxy, 'String');
minLength = get(handles.edit3_minLength, 'String');


% --------------- Check User Input ----------------

% min
if isempty( min )
    errordlg('Parameter "Minimum Frame Number" is requied by the algorithm.','Error','modal')
    return

elseif isnan(str2double(min)) || str2double(min) < 0
    errordlg('Please provide a valid value to parameter "Minimum Frame Number".','Error','modal')
    return
else
    min = str2double(min);
end    
    
% max
if isempty( max )
    errordlg('Parameter "Maximum Frame Number" is requied by the algorithm.','Error','modal')
    return

elseif isnan(str2double(max)) || str2double(max) < 0 
    errordlg('Please provide a valid value to parameter "Maximum Frame Number".','Error','modal')
    return
    
elseif str2double(max) > userData.MD.nFrames_
    errordlg('Parameter "Maximum Frame Number" can not be larger than the total number of frames.','Error','modal')
    return
    
elseif str2double(max) < min
    errordlg('Large frame number should be larger than small frame number.','Error','modal')
    return
        
else
    max = str2double(max);
end    

% filterSigma
if isempty( filterSigma )
    errordlg('Parameter "Filter Sigma" is requied by the algorithm.','Error','modal')
    return

elseif isnan(str2double(filterSigma)) || str2double(filterSigma) < 0
    errordlg('Please provide a valid value to parameter "Filter Sigma".','Error','modal')
    return
else
    filterSigma = str2double(filterSigma);
end

% dragtailLength
if isempty( dragtailLength )
    errordlg('Parameter "Dragtail Length" is requied by the algorithm.','Error','modal')
    return

elseif isnan(str2double(dragtailLength)) || str2double(dragtailLength) < 0
    errordlg('Please provide a valid value to parameter "Dragtail Length".','Error','modal')
    return
else
    dragtailLength = str2double(dragtailLength);
end  

if all(cellfun(@(x)isempty(x), {minx maxx miny maxy}))
    
    imageRange = [];
    
else

% minx
if isnan(str2double(minx)) || str2double(minx) < 0
    errordlg('Please provide a valid value to image range "min X".','Error','modal')
    return
else
    minx = (str2double(minx));
end    
    
% maxx
if isnan(str2double(maxx)) || str2double(maxx) < 0 
    errordlg('Please provide a valid value to image range "max X".','Error','modal')
    return
    
elseif str2double(maxx) > userData.MD.imSize_(1)
    errordlg('Image range "max X" can not be larger than the width of image.','Error','modal')
    return
    
elseif str2double(maxx) < minx
    errordlg('Image range "min X" cannot be larger than "max X".','Error','modal')
    return
        
else
    maxx = (str2double(maxx));
end  

% miny
if isnan(str2double(miny)) || str2double(miny) < 0
    errordlg('Please provide a valid value to image range "min Y".','Error','modal')
    return
else
    miny = (str2double(miny));
end    
    
% maxy
if isnan(str2double(maxy)) || str2double(maxy) < 0 
    errordlg('Please provide a valid value to image range "max Y".','Error','modal')
    return
    
elseif str2double(maxy) > userData.MD.imSize_(2)
    errordlg('Image range "max Y" can not be larger than the width of image.','Error','modal')
    return
    
elseif str2double(maxy) < miny
    errordlg('Image range "min Y" cannot be larger than "max Y".','Error','modal')
    return
        
else
    maxy = (str2double(maxy));
end  


imageRange = [minx maxx ; miny maxy];

end

% minLength
if isempty( minLength )
    errordlg('Parameter "Ignore tracks shorter than" is required by the algorithm.','Error','modal')
    return

elseif isnan(str2double(minLength)) || str2double(minLength) < 0
    errordlg('Please provide a valid value to parameter "Ignore tracks shorter than"".','Error','modal')
    return
else
    minLength = str2double(minLength);
end    

% -------- Set parameter --------

visualParams = userData.crtProc.visualParams_;

visualParams.otmn.startend = [min max];
visualParams.otmn.dragtailLength = dragtailLength;
visualParams.otmn.saveMovie = get(handles.checkbox3_save, 'Value');
visualParams.otmn.minLength = minLength;

if visualParams.otmn.saveMovie
    
    if isempty(filename)
        
        errordlg('Please specify a file name for result movie.','Error','modal')
        return        
    elseif isempty(path)
        
        errordlg('Please specify a path for result movie.','Error','modal')
        return
    end
    
    visualParams.otmn.movieName = filename;
    visualParams.otmn.dir2saveMovie = path;
end

visualParams.otmn.filterSigma = filterSigma;

visualParams.otmn.onlyTracks = get(handles.checkbox3_onlyTracks, 'Value');
visualParams.otmn.colorTracks = get(handles.checkbox3_colorTracks, 'Value');
visualParams.otmn.classifyGaps = get(handles.checkbox3_classifyGaps, 'Value');
visualParams.otmn.highlightES = get(handles.checkbox3_highlightES, 'Value');
visualParams.otmn.classifyLft = get(handles.checkbox3_classifyLft, 'Value');

visualParams.otmn.imageRange = imageRange;

if get(handles.radiobutton3_display_0, 'Value')
    visualParams.otmn.showRaw = 0;
    
elseif get(handles.radiobutton3_display_1, 'Value')
    visualParams.otmn.showRaw = 1;
    
elseif get(handles.radiobutton3_display_2, 'Value')
    visualParams.otmn.showRaw = 2;
end

if get(handles.radiobutton3_scale_0, 'Value')
    visualParams.otmn.intensityScale = 0;
    
elseif get(handles.radiobutton3_scale_1, 'Value')
    visualParams.otmn.intensityScale = 1;
    
elseif get(handles.radiobutton3_scale_2, 'Value')
    visualParams.otmn.intensityScale = 2;
end

userData.crtProc.setVisualParams(visualParams)
set(handles.figure1, 'UserData', userData);



overlayTracksMovieNew(userData.tracksFinal, ...
    visualParams.otmn.startend, visualParams.otmn.dragtailLength, visualParams.otmn.saveMovie, ...
    visualParams.otmn.movieName, visualParams.otmn.filterSigma, visualParams.otmn.classifyGaps, ...
    visualParams.otmn.highlightES, visualParams.otmn.showRaw, visualParams.otmn.imageRange, ...
    visualParams.otmn.onlyTracks, visualParams.otmn.classifyLft, visualParams.otmn.diffAnalysisRes, ...
    visualParams.otmn.intensityScale, visualParams.otmn.colorTracks, visualParams.otmn.firstImageFile, ...
    visualParams.otmn.dir2saveMovie,visualParams.otmn.minLength)



function edit3_min_Callback(hObject, eventdata, handles)
% hObject    handle to edit3_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3_min as text
%        str2double(get(hObject,'String')) returns contents of edit3_min as a double


% --- Executes during object creation, after setting all properties.
function edit3_min_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit3_max_Callback(hObject, eventdata, handles)
% hObject    handle to edit3_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3_max as text
%        str2double(get(hObject,'String')) returns contents of edit3_max as a double


% --- Executes during object creation, after setting all properties.
function edit3_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox3_save.
function checkbox3_save_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox3_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox3_save

userData = get(handles.figure1, 'UserData');

if get(hObject, 'Value')

    set(handles.text3_filename, 'Enable', 'on')
    set(handles.text3_mov, 'Enable', 'on')
    set(handles.pushbutton3_path, 'Enable', 'on')
    
    if isempty(get(handles.edit3_path, 'String'))
        set(handles.edit3_path, 'Enable', 'on', 'String', userData.crtProc.visualParams_.otmn.dir2saveMovie)
    else
        set(handles.edit3_path, 'Enable', 'on')
    end
    
    if isempty(get(handles.edit3_filename, 'String'))
        set(handles.edit3_filename, 'Enable', 'on', 'String', userData.crtProc.visualParams_.otmn.movieName)
    else
        set(handles.edit3_filename, 'Enable', 'on')
    end
    
else
    set(handles.text3_filename, 'Enable', 'off')
    set(handles.edit3_filename, 'Enable', 'off')
    set(handles.text3_mov, 'Enable', 'off')
    set(handles.pushbutton3_path, 'Enable', 'off')
    set(handles.edit3_path, 'Enable', 'off')    
    
end



function edit3_filename_Callback(hObject, eventdata, handles)
% hObject    handle to edit3_filename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3_filename as text
%        str2double(get(hObject,'String')) returns contents of edit3_filename as a double


% --- Executes during object creation, after setting all properties.
function edit3_filename_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3_filename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit3_filterSigma_Callback(hObject, eventdata, handles)
% hObject    handle to edit3_filterSigma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3_filterSigma as text
%        str2double(get(hObject,'String')) returns contents of edit3_filterSigma as a double


% --- Executes during object creation, after setting all properties.
function edit3_filterSigma_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3_filterSigma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit3_dragtailLength_Callback(hObject, eventdata, handles)
% hObject    handle to edit3_dragtailLength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3_dragtailLength as text
%        str2double(get(hObject,'String')) returns contents of edit3_dragtailLength as a double


% --- Executes during object creation, after setting all properties.
function edit3_dragtailLength_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3_dragtailLength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox3_highlightES.
function checkbox3_highlightES_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox3_highlightES (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox3_highlightES


% --- Executes on button press in checkbox3_onlyTracks.
function checkbox3_onlyTracks_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox3_onlyTracks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox3_onlyTracks

if get(hObject, 'Value')
   
    set(handles.checkbox3_classifyGaps, 'Value', 0, 'Enable', 'off')
    set(handles.checkbox3_highlightES, 'Value', 0, 'Enable', 'off')
    set(handles.checkbox3_classifyLft, 'Value', 0, 'Enable', 'off')
    set(handles.text44, 'Enable', 'off')
    
elseif ~get(handles.checkbox3_colorTracks, 'Value')
    
    set(handles.checkbox3_classifyGaps, 'Enable', 'on')
    set(handles.checkbox3_highlightES, 'Enable', 'on')
    set(handles.checkbox3_classifyLft, 'Enable', 'on')    
    set(handles.text44, 'Enable', 'on')
end


% --- Executes on button press in checkbox3_colorTracks.
function checkbox3_colorTracks_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox3_colorTracks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox3_colorTracks

if get(hObject, 'Value')
   
    set(handles.checkbox3_classifyGaps, 'Value', 0, 'Enable', 'off')
    set(handles.checkbox3_highlightES, 'Value', 0, 'Enable', 'off')
    set(handles.checkbox3_classifyLft, 'Value', 0, 'Enable', 'off')
    set(handles.text44, 'Enable', 'off')

elseif ~get(handles.checkbox3_onlyTracks, 'Value')
    
    set(handles.checkbox3_classifyGaps, 'Enable', 'on')
    set(handles.checkbox3_highlightES, 'Enable', 'on')
    set(handles.checkbox3_classifyLft, 'Enable', 'on')    
    set(handles.text44, 'Enable', 'on')

end


% --- Executes on button press in checkbox2_plotX.
function checkbox2_plotX_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox2_plotX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox2_plotX


% --- Executes on button press in checkbox2_plotY.
function checkbox2_plotY_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox2_plotY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox2_plotY


% --- Executes on button press in checkbox2_plotA.
function checkbox2_plotA_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox2_plotA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox2_plotA


% --- Executes on button press in checkbox2_inOneFigure.
function checkbox2_inOneFigure_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox2_inOneFigure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox2_inOneFigure


% --- Executes on button press in checkbox13.
function checkbox13_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox13


% --- Executes on button press in pushbutton2_display.
function pushbutton2_display_Callback(hObject, eventdata, handles)

% Tool 2: plotCompTrack

userData = get(handles.figure1, 'UserData');

trackid = get(handles.edit2_trackid, 'String');

% max
if isempty( trackid )
    errordlg('Parameter "Track Number" is requied by the algorithm.','Error','modal')
    return

elseif isnan(str2double(trackid)) || str2double(trackid) < 0 
    errordlg('Please provide a valid value to parameter "Track Number".','Error','modal')
    return
    
elseif str2double(trackid) > length(userData.tracksFinal)
    errordlg('Parameter "Track Number" can not be larger than the total number of tracks.','Error','modal')
    return
   
        
else
    trackid = str2double(trackid);
end   

% -------- Set parameter --------
% Tool 2: plotCompTrack

visualParams = userData.crtProc.visualParams_;

visualParams.pct.trackid = trackid;
visualParams.pct.plotX = get(handles.checkbox2_plotX, 'Value');
visualParams.pct.plotY = get(handles.checkbox2_plotY, 'Value');
visualParams.pct.plotA = get(handles.checkbox2_plotA, 'Value');
visualParams.pct.inOneFigure = get(handles.checkbox2_inOneFigure, 'Value');

userData.crtProc.setVisualParams(visualParams)
set(handles.figure1, 'UserData', userData);

plotCompTrack(userData.tracksFinal(visualParams.pct.trackid), ...
    visualParams.pct.plotX, visualParams.pct.plotY, visualParams.pct.plotA, ...
    visualParams.pct.inOneFigure, visualParams.pct.plotAggregState)




% --- Executes on button press in pushbutton3_path.
function pushbutton3_path_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

path = uigetdir(get(handles.edit3_path, 'String'), 'Select a Path...');
if path == 0
    return;
end

set(handles.edit3_path, 'String', path)

function edit3_path_Callback(hObject, eventdata, handles)
% hObject    handle to edit3_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3_path as text
%        str2double(get(hObject,'String')) returns contents of edit3_path as a double


% --- Executes during object creation, after setting all properties.
function edit3_path_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit1_min_Callback(hObject, eventdata, handles)
% hObject    handle to edit1_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1_min as text
%        str2double(get(hObject,'String')) returns contents of edit1_min as a double


% --- Executes during object creation, after setting all properties.
function edit1_min_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit1_max_Callback(hObject, eventdata, handles)
% hObject    handle to edit1_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1_max as text
%        str2double(get(hObject,'String')) returns contents of edit1_max as a double


% --- Executes during object creation, after setting all properties.
function edit1_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu1_color.
function popupmenu1_color_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1_color (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1_color contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1_color

ud = cellstr(get(hObject,'UserData'));

if strcmp(ud{get(hObject,'Value')}, '1')
    
     set(handles.text1_marker, 'Enable', 'off')
     set(handles.popupmenu1_marker, 'Enable', 'off', 'Value', length(get(handles.popupmenu1_marker, 'String')))
else
    set(handles.text1_marker, 'Enable', 'on')
    set(handles.popupmenu1_marker, 'Enable', 'on')
end

% --- Executes during object creation, after setting all properties.
function popupmenu1_color_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1_color (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox1_indicateSE.
function checkbox1_indicateSE_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1_indicateSE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1_indicateSE




% --- Executes on button press in checkbox1_flipXY.
function checkbox1_flipXY_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1_flipXY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1_flipXY


% --- Executes on button press in checkbox1_ask4sel.
function checkbox1_ask4sel_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1_ask4sel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1_ask4sel



function edit1_offset2_Callback(hObject, eventdata, handles)
% hObject    handle to edit1_offset2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1_offset2 as text
%        str2double(get(hObject,'String')) returns contents of edit1_offset2 as a double


% --- Executes during object creation, after setting all properties.
function edit1_offset2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1_offset2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit1_offset1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1_offset1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1_offset1 as text
%        str2double(get(hObject,'String')) returns contents of edit1_offset1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_offset1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1_offset1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton1_display.
function pushbutton1_display_Callback(hObject, eventdata, handles)
% Tool 1: plotTrakcs2D

userData = get(handles.figure1, 'UserData');

% --------------- Check User Input ----------------

min = get(handles.edit1_min, 'String');
max = get(handles.edit1_max, 'String');
dx = get(handles.edit1_offset1, 'String');
dy = get(handles.edit1_offset2, 'String');
filename = get(handles.edit1_image, 'String');
minLength = get(handles.edit1_minLength, 'String');

% min
if isempty( min )
    errordlg('Parameter "Minimum Frame Number" is requied by the algorithm.','Error','modal')
    return

elseif isnan(str2double(min)) || str2double(min) < 0
    errordlg('Please provide a valid value to parameter "Minimum Frame Number".','Error','modal')
    return
else
    min = str2double(min);
end    
    
% max
if isempty( max )
    errordlg('Parameter "Maximum Frame Number" is requied by the algorithm.','Error','modal')
    return

elseif isnan(str2double(max)) || str2double(max) < 0 
    errordlg('Please provide a valid value to parameter "Maximum Frame Number".','Error','modal')
    return
    
elseif str2double(max) > userData.MD.nFrames_
    errordlg('Parameter "Maximum Frame Number" can not be larger than the total number of frames.','Error','modal')
    return
    
elseif str2double(max) < min
    errordlg('Large frame number should be larger than small frame number.','Error','modal')
    return
        
else
    max = str2double(max);
end    

% dx
if isempty( dx )
    errordlg('Parameter "dX" is requied by the algorithm.','Error','modal')
    return

elseif isnan(str2double(dx)) || str2double(dx) < 0
    errordlg('Please provide a valid value to parameter "dX".','Error','modal')
    return
else
    dx = str2double(dx);
end    

% dy
if isempty( dy )
    errordlg('Parameter "dY" is requied by the algorithm.','Error','modal')
    return

elseif isnan(str2double(dy)) || str2double(dy) < 0
    errordlg('Please provide a valid value to parameter "dY".','Error','modal')
    return
else
    dy = str2double(dy);
end    

if ~isempty(filename) && ~exist(filename, 'file')
    errordlg('The image file you specified does not exist.','Error','modal')
    return    
end

% minLength
if isempty( minLength )
    errordlg('Parameter "Ignore tracks shorter than" is required by the algorithm.','Error','modal')
    return

elseif isnan(str2double(minLength)) || str2double(minLength) < 0
    errordlg('Please provide a valid value to parameter "Ignore tracks shorter than"".','Error','modal')
    return
else
    minLength = str2double(minLength);
end    

% -------- Set parameter --------

% Tool 1: plotTrakcs2D

visualParams = userData.crtProc.visualParams_;
colorUserData = get(handles.popupmenu1_color, 'UserData');
markerUserData = get(handles.popupmenu1_marker, 'UserData');

visualParams.pt2D.timeRange = [min max];
visualParams.pt2D.colorTime = colorUserData{get(handles.popupmenu1_color, 'Value')};
visualParams.pt2D.markerType = markerUserData{get(handles.popupmenu1_marker, 'Value')};

visualParams.pt2D.indicateSE = get(handles.checkbox1_indicateSE, 'Value');
visualParams.pt2D.flipXY = get(handles.checkbox1_flipXY, 'Value');
visualParams.pt2D.ask4sel = get(handles.checkbox1_ask4sel, 'Value');
visualParams.pt2D.minLength = minLength;

visualParams.pt2D.offset = [dx dy];

if isempty(filename)
   
    visualParams.pt2D.imageDir = [];
    visualParams.pt2D.image = [];    
else
    visualParams.pt2D.imageDir = filename;
    visualParams.pt2D.image = imread(filename);    
end

userData.crtProc.setVisualParams(visualParams)
set(handles.figure1, 'UserData', userData);


% visualParams.pt2D.timeRange, visualParams.pt2D.colorTime, ... % Commentable
%         visualParams.pt2D.markerType, visualParams.pt2D.indicateSE, visualParams.pt2D.newFigure,  ...
%         visualParams.pt2D.flipXY, visualParams.pt2D.ask4sel, visualParams.pt2D.offset, visualParams.pt2D.imageDir
    
% Call function
plotTracks2D(userData.tracksFinal, ...
        visualParams.pt2D.timeRange, visualParams.pt2D.colorTime, ...
        visualParams.pt2D.markerType, visualParams.pt2D.indicateSE, ...
        visualParams.pt2D.newFigure, visualParams.pt2D.image, ...
        visualParams.pt2D.flipXY, visualParams.pt2D.ask4sel,...
        visualParams.pt2D.offset,visualParams.pt2D.minLength);



% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2


% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function edit2_trackid_Callback(hObject, eventdata, handles)
% hObject    handle to edit2_trackid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2_trackid as text
%        str2double(get(hObject,'String')) returns contents of edit2_trackid as a double


% --- Executes during object creation, after setting all properties.
function edit2_trackid_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2_trackid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu1_marker.
function popupmenu1_marker_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1_marker (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1_marker contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1_marker


% --- Executes during object creation, after setting all properties.
function popupmenu1_marker_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1_marker (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton1_image.
function pushbutton1_image_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1_image (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

userData = get(handles.figure1, 'UserData');

[file,path] = uigetfile('*.*','Select an Image',...
             get(handles.edit1_image, 'String'));
        
if ~any([file,path])
    return;
end

set(handles.edit1_image, 'String', [path file])



function edit1_image_Callback(hObject, eventdata, handles)
% hObject    handle to edit1_image (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1_image as text
%        str2double(get(hObject,'String')) returns contents of edit1_image as a double


% --- Executes during object creation, after setting all properties.
function edit1_image_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1_image (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox3_classifyGaps.
function checkbox3_classifyGaps_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox3_classifyGaps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox3_classifyGaps


% --- Executes on button press in checkbox3_classifyLft.
function checkbox3_classifyLft_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox3_classifyLft (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox3_classifyLft


% --- Executes on button press in pushbutton3_crop.
function pushbutton3_crop_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3_crop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

userData = get(handles.figure1, 'UserData');

imageFileName = userData.crtProc.visualParams_.pt2D.imageDir;

% if no background image
if  isempty( imageFileName )
    imageFileName = userData.crtProc.visualParams_.otmn.firstImageFile;
end

try
    I = imread(imageFileName);   
    
catch ME
    errordlg(['Fail to open the image for cropping. Image File: ' imageFileName],'Error','modal')
    return
end

hFigure = figure('WindowStyle', 'modal', 'Name', 'Crop Image');
imshow(I)
title('Double click cropped area when finished.')
[I2, rect] = imcrop(hFigure);

if ishandle(hFigure)
    delete(hFigure)
end

if isempty(rect)
   return 
end

minY = floor(rect(1));
maxY = floor(rect(1) + rect(3));
minX = floor(rect(2));
maxX = floor(rect(2) + rect(4));

set(handles.edit3_minx, 'String', num2str(minX))
set(handles.edit3_maxx, 'String', num2str(maxX))
set(handles.edit3_miny, 'String', num2str(minY))
set(handles.edit3_maxy, 'String', num2str(maxY))




function edit3_minx_Callback(hObject, eventdata, handles)
% hObject    handle to edit3_minx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3_minx as text
%        str2double(get(hObject,'String')) returns contents of edit3_minx as a double


% --- Executes during object creation, after setting all properties.
function edit3_minx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3_minx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit3_maxx_Callback(hObject, eventdata, handles)
% hObject    handle to edit3_maxx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3_maxx as text
%        str2double(get(hObject,'String')) returns contents of edit3_maxx as a double


% --- Executes during object creation, after setting all properties.
function edit3_maxx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3_maxx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit3_miny_Callback(hObject, eventdata, handles)
% hObject    handle to edit3_miny (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3_miny as text
%        str2double(get(hObject,'String')) returns contents of edit3_miny as a double


% --- Executes during object creation, after setting all properties.
function edit3_miny_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3_miny (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit3_maxy_Callback(hObject, eventdata, handles)
% hObject    handle to edit3_maxy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3_maxy as text
%        str2double(get(hObject,'String')) returns contents of edit3_maxy as a double


% --- Executes during object creation, after setting all properties.
function edit3_maxy_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3_maxy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton3_clear.
function pushbutton3_clear_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3_clear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.edit3_minx, 'String', '');
set(handles.edit3_miny, 'String', '');
set(handles.edit3_maxx, 'String', '');
set(handles.edit3_maxy, 'String', '');



function edit3_minLength_Callback(hObject, eventdata, handles)
% hObject    handle to edit3_minLength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3_minLength as text
%        str2double(get(hObject,'String')) returns contents of edit3_minLength as a double



function edit1_minLength_Callback(hObject, eventdata, handles)
% hObject    handle to edit1_minLength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1_minLength as text
%        str2double(get(hObject,'String')) returns contents of edit1_minLength as a double
