function varargout = newMovieDataGUI(varargin)
% NEWMOVIEDATAGUI M-file for newMovieDataGUI.fig
%      NEWMOVIEDATAGUI, by itself, creates a new NEWMOVIEDATAGUI or raises the existing
%      singleton*.
%
%      H = NEWMOVIEDATAGUI returns the handle to a new NEWMOVIEDATAGUI or the handle to
%      the existing singleton*.
%
%      NEWMOVIEDATAGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in NEWMOVIEDATAGUI.M with the given input arguments.
%
%      NEWMOVIEDATAGUI('Property','Value',...) creates a new NEWMOVIEDATAGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before newMovieDataGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to newMovieDataGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help newMovieDataGUI

% Last Modified by GUIDE v2.5 17-Nov-2010 14:29:25

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @newMovieDataGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @newMovieDataGUI_OutputFcn, ...
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


% --- Executes just before newMovieDataGUI is made visible.
function newMovieDataGUI_OpeningFcn(hObject, eventdata, handles, varargin)
%
% newMovieDataGUI('mainFig', handles.figure1) - call from movieSelector
% newMovieDataGUI('mainFig', handles.figure1, 'overview', MD) - call movie data
% viewer from movieSelecter 'Detail ...'
%
% Useful tools:
%
% User Data:
% 
% userData.channels - array of Channel objects
% userData.mainFig - handle of movie selector GUI
% userData.handles_main - 'handles' of movie selector GUI
%
% userData.setChannelFig - handle of channel set-up figure
% userData.iconHelpFig - handle of help dialog
%
% NOTE: If newMovieDataGUI is under the "Overview" mode, additionally, 
% 
% userData.MD - the handle of selected MovieData object
%
%

[copyright openHelpFile] = userfcn_softwareConfig(handles);
set(handles.text_copyright, 'String', copyright)

userData = get(handles.figure1, 'UserData');
% Choose default command line output for newMovieDataGUI
handles.output = hObject;

% Set channel object array
userData.channels = [];

% Load help icon from dialogicons.mat
load lccbGuiIcons.mat
supermap(1,:) = get(hObject,'color');

userData.colormap = supermap;
userData.questIconData = questIconData;

axes(handles.axes_help);
Img = image(questIconData);
set(hObject,'colormap',supermap);
set(gca, 'XLim',get(Img,'XData'),'YLim',get(Img,'YData'),...
    'visible','off');
set(Img,'ButtonDownFcn',@icon_ButtonDownFcn);

if openHelpFile
    set(Img, 'UserData', struct('class', 'NewMovieData'))
end

if nargin > 3
    
    t = find(strcmp(varargin,'mainFig'));
    if isempty(t)
        error('User-defined: input error, correct statement: setupMovieDataGUI(''mainFig'',handles.figure1)');
    end
    userData.mainFig = varargin{t+1};
    userData.handles_main = guidata(userData.mainFig);
    
    t = find(strcmp(varargin, 'overview'));
    
    if ~isempty(t)
        % The panel is launched as a movie data viewer

        % Set userData.MD - MovieData object from movieSelectorGUI
        assert( isa(varargin{t(1)+1}, 'MovieData'), 'User-defined: invalid movie data input.' )

        userData.MD = varargin{t(1)+1};
        userData.channels = userData.MD.channels_;
        
        % Channel listbox
        cPath = cell(1,length(userData.channels));
        for i = 1: length(userData.channels)
            cPath{i} = userData.channels(i).channelPath_;
        end
        set(handles.listbox_channel, 'String', cPath)
        
        % GUI setting
        set(handles.pushbutton_delete, 'Visible', 'off')
        set(handles.pushbutton_add, 'Visible', 'off')
        set(handles.pushbutton_cancel, 'Visible', 'off')
        set(handles.pushbutton_done, 'String', 'Ok')
        set(handles.pushbutton_setting_chan, 'String', 'Advanced Channel Setting (Editable)')
        set(handles.pushbutton_output, 'Visible', 'off')
        
        set(hObject, 'Name', 'Movie Detail')
        set(handles.edit_output, 'Enable', 'on')
        set(handles.uipanel_notes, 'Title', 'Notes (Editable and will be saved)')
        set(handles.text_body_1, 'Visible', 'on')
        set(handles.edit_path, 'Visible','on', 'String', [userData.MD.movieDataPath_ userData.MD.movieDataFileName_])
        
        % GUI setting - parameters
        set(handles.edit_output, 'String', userData.MD.outputDirectory_)
        set(handles.edit_ps, 'Enable', 'inactive', 'String', userData.MD.pixelSize_)
        set(handles.edit_ti, 'Enable', 'inactive', 'String', userData.MD.timeInterval_)
        set(handles.edit_na, 'Enable', 'inactive', 'String', userData.MD.numAperature_)
        set(handles.edit_cb, 'Enable', 'inactive', 'String', userData.MD.camBitdepth_)
        
        set(handles.edit_notes, 'String', userData.MD.notes_)
    end
    
    
% else
%     error('User-defined: For now this new movie data GUI can only be called from movieDataSelector.')
    
end


% Update handles structure
set(handles.figure1,'UserData',userData)
guidata(hObject, handles);

% UIWAIT makes newMovieDataGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = newMovieDataGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton_cancel.
function pushbutton_cancel_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
delete(handles.figure1);


% --- Executes on button press in pushbutton_done.
function pushbutton_done_Callback(hObject, eventdata, handles)

userData = get(handles.figure1,'UserData');
userData_main = get(userData.mainFig, 'UserData');

% If newMovieDataGUI is under "Overview" mode
notes = get(handles.edit_notes, 'string');

if isfield(userData, 'MD')
    
    userData.MD.setNotes(notes)
    delete(handles.figure1)
    return
end

% Verify channels are given
if ~isfield(userData, 'channels') || isempty(userData.channels)
    errordlg('Please provide at least one channel path.',...
        'Empty Channel','modal');
    return;    
end

if ~isa(userData.channels(1), 'Channel')
   error('User-defined: userData.channels are not of class ''Channel''') 
end

% Check output path
outputDir = get(handles.edit_output, 'String');

if isempty(outputDir)
    errordlg('Please provide an output path to save your results.', ...
               'Empty Output Path', 'modal');
    return;    
end

if ~exist(outputDir, 'dir')
    errordlg('Please provide an valid path to save your results.', ...
               'Invalid Output Path', 'modal');
    return; 
end


% Assign user-defined parameters

pixelSize = get(handles.edit_ps, 'string');
timeInterval = get(handles.edit_ti, 'string');
numAperature = get(handles.edit_na, 'string');
camBitdepth = get(handles.edit_cb, 'string');

% Creat the movieData object (try/catch block) -------------!!!
try
    MD = MovieData(userData.channels, outputDir, [], [],...
        notes, pixelSize, timeInterval, numAperature, camBitdepth);
catch ME
    errormsg = sprintf([ME.message '.\n\nCreating movie data failed.']);
    errordlg(errormsg, 'User Input Error','modal');
    return;
end  

try
%Sanity Check
    MD.sanityCheck

catch ME
    delete(MD);
    errormsg = sprintf('%s.\n\nPlease check your movie data. Movie data is not saved.',ME.message);
    errordlg(errormsg,'Channel Error','modal');
    return;
end

% Ask user where to save the movie data file
[file,path] = uiputfile('*.mat','Find a place to save your movie data',...
             [outputDir filesep 'movieData.mat']);
        
if ~any([file,path])
    return;
end

% Check if files in movie list are saved in the same file
contentlist = get(userData.handles_main.listbox_movie, 'String');

if any(strcmp([path file], contentlist))
    errordlg('Cannot overwrite a movie data file which is already in the movie list. Please choose another file name or another path.','Error','modal');
    return
end

% After checking file directory, set directory to movie data
MD.setMovieDataPath(path)
MD.setMovieDataFileName(file)

% Save Movie Data to disk
MD.saveMovieData

% Save Movie Data to movie selector panel
userData_main.MD = cat(2, userData_main.MD, MD);

% Refresh movie list box in movie selector panel
contentlist{end+1} = [path file];
set(userData.handles_main.listbox_movie, 'String', contentlist, 'Value', length(contentlist))
title = sprintf('Movie List: %s/%s movie(s)', num2str(length(contentlist)), num2str(length(contentlist)));
set(userData.handles_main.text_movie_1, 'String', title)

set(userData.mainFig, 'UserData', userData_main)

% Delete current window
delete(handles.figure1)




function edit_ti_Callback(hObject, eventdata, handles)
% hObject    handle to edit_ti (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_ti as text
%        str2double(get(hObject,'String')) returns contents of edit_ti as a double


% --- Executes during object creation, after setting all properties.
function edit_ti_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_ti (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_ps_Callback(hObject, eventdata, handles)
% hObject    handle to edit_ps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_ps as text
%        str2double(get(hObject,'String')) returns contents of edit_ps as a double


% --- Executes during object creation, after setting all properties.
function edit_ps_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_ps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox_channel.
function listbox_channel_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_channel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_channel contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_channel


% --- Executes during object creation, after setting all properties.
function listbox_channel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_channel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox_channel controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_delete.
function pushbutton_delete_Callback(hObject, eventdata, handles)

userData = get(handles.figure1, 'Userdata');

contents = get(handles.listbox_channel,'String');
% Return if list is empty
if isempty(contents)
    return;
end
num = get(handles.listbox_channel,'Value');

% Delete channel object
delete(userData.channels(num))
userData.channels(num) = [];

% Refresh listbox_channel
contents(num) = [ ];
set(handles.listbox_channel,'String',contents);

% Point 'Value' to the second last item in the list once the 
% last item has been deleted
if num>length(contents) && num>1
    set(handles.listbox_channel,'Value',length(contents));
end

set(handles.figure1, 'Userdata', userData)
guidata(hObject, handles);

% --- Executes on button press in pushbutton_add.
function pushbutton_add_Callback(hObject, eventdata, handles)

set(handles.listbox_channel, 'Value', 1)

userData = get(handles.figure1, 'UserData');
userData_main = get(userData.handles_main.figure1, 'UserData');

path = uigetdir(userData_main.userDir, 'Add Channels ...');
if path == 0
    return;
end
% Input validation function ...
% Get current list
contents = get(handles.listbox_channel,'String');
if any(strcmp(contents,path))
   warndlg('This directory has been selected! Please select a differenct directory.',...
       'Warning','modal');
   return; 
end

% Create path object and save it to userData
userData.channels = cat(2, userData.channels, Channel(path));

% Refresh listbox_channel
contents{end+1} = path;
set(handles.listbox_channel,'string',contents);

% Set user directory
sepDir = regexp(path, filesep, 'split');
dir = sepDir{1};
for i = 2: length(sepDir)-1
    dir = [dir filesep sepDir{i}];
end
userData_main.userDir = dir;

set(handles.figure1, 'Userdata', userData)
set(userData.handles_main.figure1, 'UserData', userData_main)
guidata(hObject, handles);


% --- Executes during object deletion, before destroying properties.
function figure1_DeleteFcn(hObject, eventdata, handles)

userData = get(handles.figure1, 'UserData');

if isfield(userData, 'iconHelpFig') && ishandle(userData.iconHelpFig)
   delete(userData.iconHelpFig) 
end




function edit_notes_Callback(hObject, eventdata, handles)
% hObject    handle to edit_notes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_notes as text
%        str2double(get(hObject,'String')) returns contents of edit_notes as a double


% --- Executes during object creation, after setting all properties.
function edit_notes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_notes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_output_Callback(hObject, eventdata, handles)
% hObject    handle to edit_output (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_output as text
%        str2double(get(hObject,'String')) returns contents of edit_output as a double


% --- Executes during object creation, after setting all properties.
function edit_output_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_output (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_output.
function pushbutton_output_Callback(hObject, eventdata, handles)

pathname = uigetdir(pwd);
if isnumeric(pathname)
    return;
end

set(handles.edit_output, 'String', pathname);


% --- Executes on button press in pushbutton_setting_chan.
function pushbutton_setting_chan_Callback(hObject, eventdata, handles)

userData = get(handles.figure1, 'UserData');

if isempty(userData.channels) 
    return
end
assert( isa(userData.channels(1), 'Channel'), 'User-defined: Not a valid ''Channel'' object');

userData.setChannelFig = newChannelGUI('mainFig', handles.figure1, 'modal');

set(handles.figure1,'UserData',userData);


function edit_na_Callback(hObject, eventdata, handles)
% hObject    handle to edit_na (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_na as text
%        str2double(get(hObject,'String')) returns contents of edit_na as a double


% --- Executes during object creation, after setting all properties.
function edit_na_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_na (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_cb_Callback(hObject, eventdata, handles)
% hObject    handle to edit_cb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_cb as text
%        str2double(get(hObject,'String')) returns contents of edit_cb as a double


% --- Executes during object creation, after setting all properties.
function edit_cb_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_cb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



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
