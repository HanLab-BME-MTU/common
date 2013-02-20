function varargout = detectionVisualGUI(varargin)
% DETECTIONVISUALGUI M-file for detectionVisualGUI.fig
%      DETECTIONVISUALGUI, by itself, creates a new DETECTIONVISUALGUI or raises the existing
%      singleton*.
%
%      H = DETECTIONVISUALGUI returns the handle to a new DETECTIONVISUALGUI or the handle to
%      the existing singleton*.
%
%      DETECTIONVISUALGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DETECTIONVISUALGUI.M with the given input arguments.
%
%      DETECTIONVISUALGUI('Property','Value',...) creates a new DETECTIONVISUALGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before detectionVisualGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to detectionVisualGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help detectionVisualGUI

% Last Modified by GUIDE v2.5 17-Dec-2010 16:04:17

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @detectionVisualGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @detectionVisualGUI_OutputFcn, ...
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


% --- Executes just before detectionVisualGUI is made visible.
function detectionVisualGUI_OpeningFcn(hObject, eventdata, handles, varargin)

% userData.showDetFig = detectionVisualGUI('mainFig',handles.figure1, procID);
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
%       userData.firstframe - the first frame who has detection result
%       userData.lastframe - the last frame who has detection result
%       userData.file - maximum frame number
%

set(handles.text_copyright, 'String', getLCCBCopyright());

userData = get(handles.figure1, 'UserData');
% Choose default command line output for detectionVisualGUI
handles.output = hObject;

set(handles.uipanel_display, 'SelectionChangeFcn', @uipanel_display_SelectionChangeFcn)
set(handles.uipanel_scale, 'SelectionChangeFcn', @uipanel_scale_SelectionChangeFcn)

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

% Make sure output exists
chan = [];
for i = 1:length(userData.MD.channels_)
    if userData.crtProc.checkChannelOutput(i)
        chan = i; 
        break
    end
end

if isempty(chan)
   error('User-defined: the process does not have output yet.') 
end

load(userData.crtProc.outFilePaths_{chan},'movieInfo');
% Make sure detection output is valid
firstframe = [];
for i = 1:length(movieInfo)
   
    if ~isempty(movieInfo(i).amp)
        firstframe = i;
        break
    end
end

if isempty(firstframe)
   error('User-defined: there is no detection information in the output variable.') 
end


userData.movieInfo=movieInfo;
userData.chan = chan;
userData.firstframe = firstframe;
userData.lastframe = length(userData.movieInfo);

userData.toolName = {'overlayFeaturesMovie'};

% ------------------ Update input parameters ------------------------

visualParams = userData.crtProc.visualParams_;
% parameter: startend
visualParams.startend = [userData.firstframe, userData.lastframe]; 
% parameter: firstImageFile
% file = userData.MD.getImageFileNames(chan);
% visualParams.firstImageFile = [userData.MD.channels_(chan).channelPath_ filesep file{1}{1}];
% parameter: dir2saveMovie
if any(strcmp(visualParams.dir2saveMovie, {userData.MD.channels_(:).channelPath_}))
    visualParams.dir2saveMovie = userData.MD.channels_(chan).channelPath_;
end
userData.crtProc.setVisualParams(visualParams)

% -------------------- Set Parameters -------------------------------

visualParams = userData.crtProc.visualParams_;

% text_movie
str = userData.MD.getChannelPaths(chan);
limit = 60;
if length(str)>limit
    str = ['... ' str(end-limit:end)];
end
set(handles.text_movie, 'String', str)

% text_framenum
set(handles.text_framenum, 'String', ['( Available detection from frame ',num2str(userData.firstframe),' to ',num2str(userData.lastframe),' )'])

set(handles.edit_min, 'String', num2str(visualParams.startend(1)))
set(handles.edit_max, 'String', num2str(visualParams.startend(2)))
set(handles.edit_filtersigma, 'String', num2str(visualParams.filterSigma))

% Set default movie name as the name of result MAT file
% if isempty(visualParams.movieName)
[~,visualParams.movieName]=fileparts(userData.crtProc.outFilePaths_{1,chan});
userData.crtProc.setVisualParams(visualParams)
visualParams = userData.crtProc.visualParams_;
% end

if ~visualParams.saveMovie
    
    set(handles.checkbox_save, 'Value', 0)
    set(handles.text_filename, 'Enable', 'off')
    set(handles.edit_filename, 'Enable', 'off', 'String', '')
    set(handles.text_mov, 'Enable', 'off')
    set(handles.pushbutton_path, 'Enable', 'off')
    set(handles.edit_path, 'Enable', 'off')
    
else
    set(handles.edit_filename, 'String', visualParams.movieName)
    set(handles.edit_path, 'String', visualParams.dir2saveMovie)
    
end

switch visualParams.showRaw
    
    case 0 
        set(handles.radiobutton_display_0, 'Value', 1)
    case 1
        set(handles.radiobutton_display_1, 'Value', 1)
    case 2
        set(handles.radiobutton_display_2, 'Value', 1)
    otherwise
        error('User-defined: the parameter is incorrect.')
end

switch visualParams.intensityScale
    
    case 0 
        set(handles.radiobutton_scale_0, 'Value', 1)
    case 1
        set(handles.radiobutton_scale_1, 'Value', 1)
    case 2
        set(handles.radiobutton_scale_2, 'Value', 1)
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
    set(Img,'ButtonDownFcn',@icon_ButtonDownFcn,...
        'UserData', struct('class', userData.toolName{i}))

end


% Update user data and GUI data
set(hObject, 'UserData', userData);
guidata(hObject, handles);

% UIWAIT makes detectionVisualGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = detectionVisualGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1_display.
function pushbutton1_display_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1_display (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

userData = get(handles.figure1, 'UserData');


% --------------- Check User Input ----------------

min = get(handles.edit_min, 'String');
max = get(handles.edit_max, 'String');
filterSigma = get(handles.edit_filtersigma, 'String');
filename = get(handles.edit_filename, 'String');
path = get(handles.edit_path, 'String');

% min
if isempty( min )
    errordlg('Parameter "Frames to Included in Movie" is requied by the algorithm.','Error','modal')
    return

elseif isnan(str2double(min)) || str2double(min) < 0
    errordlg('Please provide a valid value to parameter "Frames to Included in Movie".','Error','modal')
    return
else
    min = str2double(min);
end    
    
% max
if isempty( max )
    errordlg('Parameter "Frames to Included in Movie" is requied by the algorithm.','Error','modal')
    return

elseif isnan(str2double(max)) || str2double(max) < 0 
    errordlg('Please provide a valid value to parameter "Frames to Included in Movie".','Error','modal')
    return
    
elseif str2double(max) > userData.MD.nFrames_
    errordlg('Parameter "Frames to Included in Movie" can not be larger than the total number of frames.','Error','modal')
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

% -------- Set parameter --------

visualParams = userData.crtProc.visualParams_;

visualParams.startend = [min max];
visualParams.saveMovie = get(handles.checkbox_save, 'Value');

if visualParams.saveMovie
    
    if isempty(filename) || isempty(path)
       errordlg('Please specify a file name and directory to save the Quick time movie.','Error','modal')
       return
    end
    
    visualParams.movieName = filename;
    visualParams.dir2saveMovie = path;
end

visualParams.filterSigma = filterSigma;

if get(handles.radiobutton_display_0, 'Value')
    visualParams.showRaw = 0;
    
elseif get(handles.radiobutton_display_1, 'Value')
    visualParams.showRaw = 1;
    
elseif get(handles.radiobutton_display_2, 'Value')
    visualParams.showRaw = 2;
end

if get(handles.radiobutton_scale_0, 'Value')
    visualParams.intensityScale = 0;
    
elseif get(handles.radiobutton_scale_1, 'Value')
    visualParams.intensityScale = 1;
    
elseif get(handles.radiobutton_scale_2, 'Value')
    visualParams.intensityScale = 2;
end
    
userData.crtProc.setVisualParams(visualParams)
set(handles.figure1, 'UserData', userData);

% Display result

overlayFeaturesMovie(userData.movieInfo, ...
    visualParams.startend, visualParams.saveMovie,visualParams.movieName,...
    visualParams.filterSigma,visualParams.showRaw,visualParams.intensityScale, ...
    visualParams.firstImageFile, visualParams.dir2saveMovie)





function edit_min_Callback(hObject, eventdata, handles)
% hObject    handle to edit_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_min as text
%        str2double(get(hObject,'String')) returns contents of edit_min as a double


% --- Executes during object creation, after setting all properties.
function edit_min_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_max_Callback(hObject, eventdata, handles)
% hObject    handle to edit_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_max as text
%        str2double(get(hObject,'String')) returns contents of edit_max as a double


% --- Executes during object creation, after setting all properties.
function edit_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_save.
function checkbox_save_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_save
userData = get(handles.figure1, 'UserData');

if get(hObject, 'Value')

    set(handles.text_filename, 'Enable', 'on')
    
    set(handles.text_mov, 'Enable', 'on')
    set(handles.pushbutton_path, 'Enable', 'on')
    if isempty(get(handles.edit_path, 'String'))
        set(handles.edit_path, 'Enable', 'on', 'String',userData.crtProc.visualParams_.dir2saveMovie)
    else
        set(handles.edit_path, 'Enable', 'on')
    end
    
    if isempty(get(handles.edit_filename, 'String'))
        set(handles.edit_filename, 'Enable', 'on', 'String', userData.crtProc.visualParams_.movieName)
    else
        set(handles.edit_filename, 'Enable', 'on')
    end
else
    set(handles.text_filename, 'Enable', 'off')
    set(handles.edit_filename, 'Enable', 'off')
    set(handles.text_mov, 'Enable', 'off')
    set(handles.pushbutton_path, 'Enable', 'off')
    set(handles.edit_path, 'Enable', 'off')    
    
end



function edit_filename_Callback(hObject, eventdata, handles)
% hObject    handle to edit_filename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_filename as text
%        str2double(get(hObject,'String')) returns contents of edit_filename as a double


% --- Executes during object creation, after setting all properties.
function edit_filename_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_filename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_filtersigma_Callback(hObject, eventdata, handles)
% hObject    handle to edit_filtersigma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_filtersigma as text
%        str2double(get(hObject,'String')) returns contents of edit_filtersigma as a double


% --- Executes during object creation, after setting all properties.
function edit_filtersigma_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_filtersigma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_path.
function pushbutton_path_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

path = uigetdir(get(handles.edit_path, 'String'), 'Select a Path...');
if path == 0
    return;
end

set(handles.edit_path, 'String', path)



function edit_path_Callback(hObject, eventdata, handles)
% hObject    handle to edit_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_path as text
%        str2double(get(hObject,'String')) returns contents of edit_path as a double


% --- Executes during object creation, after setting all properties.
function edit_path_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function uipanel_display_SelectionChangeFcn(hObject, eventdata)

function uipanel_scale_SelectionChangeFcn(hObject, eventdata)


% --- Executes on button press in pushbutton_done.
function pushbutton_done_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_done (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
delete(handles.figure1)
