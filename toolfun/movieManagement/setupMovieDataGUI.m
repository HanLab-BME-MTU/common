function varargout = setupMovieDataGUI(varargin)
% SETUPMOVIEDATAGUI M-file for SETUPMOVIEDATAGUI.fig
%      SETUPMOVIEDATAGUI, by itself, creates a new SETUPMOVIEDATAGUI or
%      raises the existing singleton*.
%
%      H = SETUPMOVIEDATAGUI returns the handle to a new SETUPMOVIEDATAGUI
%      or the handle to the existing singleton*.
%
%      SETUPMOVIEDATAGUI('CALLBACK',hObject,eventData,handles,...) calls
%      the local function named CALLBACK in SETUPMOVIEDATAGUI.M with the
%      given input arguments. 
%
%      SETUPMOVIEDATAGUI('Property','Value',...) creates a new
%      SETUPMOVIEDATAGUI or raises the existing singleton*.  Starting from
%      the left, property value pairs are applied to the GUI before
%      setupMovieDataGUI_OpeningFcn gets called.  An unrecognized property
%      name or invalid value makes property application stop.  All inputs
%      are passed to setupMovieDataGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help setupMovieDataGUI

% Last Modified by GUIDE v2.5 20-Jun-2010 21:35:25

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @setupMovieDataGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @setupMovieDataGUI_OutputFcn, ...
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


% --- Executes just before setupMovieDataGUI is made visible.
function setupMovieDataGUI_OpeningFcn(hObject, eventdata, handles, varargin) %#ok<INUSL>
%
% Useful tools:
% 
% User Data:
%
%   userData.MD - save new or loaded Movie Data object
%   userData.userDir - save the default directory when user click 'Add Channel'
%   userData.mainFig - handle of package GUI (if figure is called by package GUI)
%   userData.helpFig - handle of help dialog
%
%   userData.firstPackageName - the package class name 
%   userData.furstPackageGUI - the name of package GUI
%

userData = get(handles.figure1, 'UserData');

% Choose default command line output for setupMovieDataGUI
handles.output = hObject;



% Set callback function of radio button group uipanel_1
set(handles.uipanel_1, 'SelectionChangeFcn', @uipanel_1_SelectionChangeFcn);

% Set radio button unchecked 
set(handles.uipanel_1, 'SelectedObject', []);
 
% Set GUI data used to save MovieData object
userData.MD = [ ];

userData.userDir = pwd;

% Indicate the first package control panel the GUI will go to after 
% MovieData is sucessfully created and MovieDate setup panel is destroied.
% For future needs of multiple packages or multiple movie data input. 
% These variables can be re-defined as input variables when MovieData 
% setup panel is created.

userData.firstPackageName = 'BioSensorsPackage';
userData.firstPackageGUI = @biosensorsPackageGUI;

% Load help icon from dialogicons.mat
load lccbGuiIcons.mat
supermap(1,:) = get(hObject,'color');

userData.colormap = supermap;
userData.questIconData = questIconData;

axes(handles.axes_1);
Img = image(questIconData);
set(hObject,'colormap',supermap);
set(gca, 'XLim',get(Img,'XData'),'YLim',get(Img,'YData'),...
    'visible','off');
set(Img,'ButtonDownFcn',@icon_ButtonDownFcn);
set(Img, 'UserData', MovieData.getHelp)

% Save userdata
set(handles.figure1,'UserData',userData);

% If GUI is called by package control panel (called by menu 'New' or 'Open')
if nargin > 3
    t = find(strcmp(varargin,'mainFig'));
    if isempty(t)
        error('User-defined: input error, correct statement: setupMovieDataGUI(''mainFig'',handles.figure1)');
    end
    userData.mainFig = varargin{t+1};
    
    % Notify control panel that setupMovieData GUI is open
    setappdata(userData.mainFig, 'setupMovieDataFlag', 1);
    
    % 'Continue to package' checkbox set to invisible, which is an obvious
    % case
    set(handles.checkbox_continue, 'Visible', 'off');

end

% Save userdata
set(handles.figure1,'UserData',userData);
guidata(hObject, handles);

% UIWAIT makes setupMovieDataGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = setupMovieDataGUI_OutputFcn(hObject, eventdata, handles)  %#ok<INUSL>
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in pushbutton_done.
function pushbutton_done_Callback(hObject, eventdata, handles) %#ok<*DEFNU>
% hObject    handle to pushbutton_done (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get User Data
userData = get(handles.figure1,'UserData');

if isempty (get(handles.uipanel_1, 'SelectedObject'));
    errordlg('Please select a method to choose input channels.',...
        'Please select a method','modal');
    return;
end

% Get notes
notes = get(handles.edit_notes, 'string');

switch get(get(handles.uipanel_1, 'SelectedObject'), 'tag')
    
    case 'radiobutton_1'
        
        % Get input parameters channel path (cell array), pixel size and 
        % time interval
        channelPath = get(handles.listbox, 'string');
        pixelSize = get(handles.edit_ps, 'string');
        timeInterval = get(handles.edit_ti, 'string');
        
        
        % Check output directory 
        outputDir = get(handles.edit_output, 'String');

        if isempty(outputDir)
            errordlg('Please provide an output path to save your results.', ...
                    'User Input Error', 'modal');
            return;    
        end

        if ~exist(outputDir, 'dir')
            errordlg('Please provide an valid path to save your results.', ...
                    'User Input Error', 'modal');
            return; 
        end
        
        % Check if text box is empty
        if isempty(channelPath)
            errordlg('Please provide at least one channel path.',...
                       'User Input Error','modal');
            return;
        end
        
        if isempty(pixelSize) || ...
                isempty(timeInterval)
            errordlg('Please provide pixel size and time interval values.',...
                'User Input Error','modal');
            return;
        end
        
        % Ask the user where to save the movie data
        [file,path] = uiputfile('*.mat','Find a place to save your movie data',...
           [outputDir filesep 'movieData.mat']);
        
        if ~any([file,path])
            return;
        end
        
        % Creat the movieData object (try/catch block) -------------!!!
        try
            MD = MovieData(channelPath,pixelSize,timeInterval,path,file,notes,outputDir);
        catch ME
            errormsg = sprintf([ME.message '.\n\nMovie data is not saved.']);
            errordlg(errormsg, 'User Input Error','modal');
            return;
        end   

        try
        % Sanity Check
            MD.sanityCheck(path, file, false);
        % Catch: exception occurs - data is not saved.
        catch ME
            delete(MD);
            userData.MD = [ ];
            errormsg = sprintf('%s.\n\nPlease check your movie data. Movie data is not saved.',ME.message);
            errordlg(errormsg,'Channel Path & Movie Data Error','modal');
            return;
        end
        
        % Set MovieDate as GUI data
        userData.MD = MD;
        set(handles.figure1,'UserData',userData);
        guidata(hObject, handles);

    case 'radiobutton_2'
        
        % Movie Data sanity check has be done in callback fcn of pushbutton_open 
        % Go to the first package
        if isempty(userData.MD)
            errordlg(['Movie data has not been successfully loaded. Please' ...
                ' click ''Open File...'' button to select a valid MAT file.'],...
                'Movie Data Not Loaded','modal');
            return;
        end
        
        MD = userData.MD;
        MD.setNotes(notes)
        
    otherwise
        error('User-defined: unexpected radio button has been selected');
end


% Save MovieData
save([MD.movieDataPath_ MD.movieDataFileName_], 'MD');

% If setupMovieDataGUI is called from package panel, and package GUI has
% Movie Data loaded
if isfield(userData, 'mainFig') && ...
        isfield(get(userData.mainFig, 'Userdata'), 'MD')
    
    userData_main = get(userData.mainFig, 'Userdata');
    MD_main = userData_main.MD;
    
    user_response = questdlg(['Before loading the new movie data, do you want to save the current progress to ',MD_main.movieDataFileName_,'?'], ...
    'BioSensors Package Control Panel');

    switch lower(user_response)
        case 'yes'
            save([MD_main.movieDataPath_ MD_main.movieDataFileName_], 'MD_main');
            
            setappdata(userData.mainFig, 'setupMovieDataFlag', 0);
            delete(userData.mainFig);
            userData.firstPackageGUI(MD);
            
            delete(handles.figure1);
            
        case 'no'
            setappdata(userData.mainFig, 'setupMovieDataFlag', 0);
            delete(userData.mainFig);
            userData.firstPackageGUI(MD);    
            
            delete(handles.figure1);
            
        case 'cancel'
            % Save user data
            set(handles.figure1,'UserData', userData);
            guidata(hObject, handles);
    end
    
elseif get(handles.checkbox_continue, 'Value')

    userData.firstPackageGUI(MD);
    delete(handles.figure1)
else
    delete(handles.figure1)
end




% --- Executes on selection change in listbox.
function listbox_Callback(hObject, eventdata, handles)
% hObject    handle to listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox


% --- Executes during object creation, after setting all properties.
function listbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_del.
function pushbutton_del_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_del (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
contents = get(handles.listbox,'String');
% Return if list is empty
if isempty(contents)
    return;
end
num = get(handles.listbox,'Value');

% Delete selected item
contents(num) = [ ];

% Refresh listbox
set(handles.listbox,'String',contents);
% Point 'Value' to the second last item in the list once the 
% last item has been deleted
if (num>length(contents) && num>1)
    set(handles.listbox,'Value',length(contents));
end

guidata(hObject, handles);

% --- Executes on button press in pushbutton_add.
function pushbutton_add_Callback(hObject, eventdata, handles)

set(handles.listbox, 'Value', 1)

userData = get(handles.figure1, 'Userdata');

path = uigetdir(userData.userDir, 'Add Channels ...');
if path == 0
    return;
end
% Input validation function ...
% Get current list
contents = get(handles.listbox,'String');
if any(strcmp(contents,path))
   warndlg('This directory has been selected! Please select a differenct directory.',...
       'Warning','modal');
   return; 
end
% Add current formula to the listbox
contents{end+1} = path;
set(handles.listbox,'string',contents);

% Set user directory
sepDir = regexp(path, filesep, 'split');
dir = sepDir{1};
for i = 2: length(sepDir)-1
    dir = [dir filesep sepDir{i}];
end
userData.userDir = dir;

set(handles.figure1, 'Userdata', userData)
guidata(hObject, handles);



function edit_mat_Callback(hObject, eventdata, handles)
% hObject    handle to edit_mat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_mat as text
%        str2double(get(hObject,'String')) returns contents of edit_mat as a double


% --- Executes during object creation, after setting all properties.
function edit_mat_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_mat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_open.
function pushbutton_open_Callback(hObject, eventdata, handles)

userData = get(handles.figure1, 'UserData');
[filename, pathname] = uigetfile('*.mat','Select Movie Data MAT-file');
if ~any([filename pathname])
    return;
end
try
    pre = whos('-file', [pathname filename]);  % - Exception: fail to access .mat file
catch ME
    errordlg(ME.message,'Fail to open MAT file.','modal');
    return;
end
% Find MovieData object in .mat file before loading it
structMD = pre( logical(strcmp({pre(:).class},'MovieData')) );
switch length(structMD)
    case 0
        % Exception: No MovieData object is in the MAT file
        errordlg('No MovieData object is found in selected MAT-file.',...
            'MAT File Error','modal');
        return;
    case 1
        % Load MovieData object in .mat file
        load([pathname filename],'-mat',structMD.name);
        % Name/Rename MovieData obj as 'MD'
        eval(['MD =' structMD.name ';']);
    otherwise
        % Exception - TODO: multiple movies saved in one MAT file
        errordlg('The selected MAT-file contains more than one Movie Data.',...
            'MAT File Error','modal');
        return;
end
% Set parameter values in .mat file to GUI        
set(handles.edit_mat,'String',[pathname,filename]);
set(handles.listbox_2,'String',MD.channelPath_);
set(handles.edit_ps2,'String',MD.pixelSize_);
set(handles.edit_ti2,'String',MD.timeInterval_);
set(handles.edit_notes,'String',MD.notes_);

set(handles.edit_output,'String',MD.outputDirectory_);

try
% Sanity check of saved movie data 
MD.sanityCheck( pathname,filename,false );       
catch ME
    % If don't pass sanity check, set MovieData unloaded
    delete(MD);
    userData.MD = [ ];
    errordlg([ME.message 'Movie data is not successfully loaded.'],...
                'Channel Path & Movie Data Error','modal');   
    % Erase success notice
    set(handles.text_body8,'Visible','off');
    set(handles.figure1,'UserData',userData);
    guidata(hObject, handles);
    
    return;
end
        % Success Notice: Movie Data has been loaded
set(handles.text_body8,'visible','on');

userData.MD = MD;

set(handles.figure1,'UserData',userData);
guidata(hObject,handles);
        
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

function uipanel_1_SelectionChangeFcn(hObject, eventdata)
% Call back function of ration button group uipanel_1
handles = guidata(hObject); 

% Highlight the content under new radiobutton
switch get(eventdata.NewValue,'Tag')   % Get Tag of selected object
    case 'radiobutton_1'
      %execute this code when fontsize08_radiobutton is selected
      userfcn_fade(handles,'method1',1);
    case 'radiobutton_2'
      %execute this code when fontsize12_radiobutton is selected
      userfcn_fade(handles,'method2',1);
    otherwise
       disp('User-defined Warning: No radio button tag is ',...
           'found when SelectionChangeFcn is triggered.');
       return;
end
% Fade the content under old radiobutton
% For the first time of radiobutton selection, no need to take any action
if isempty (eventdata.OldValue)
    
    % Enable 'Finish' button
    set(handles.pushbutton_done, 'enable','on');
    return; 
end
switch get(eventdata.OldValue,'Tag')
    case 'radiobutton_1'
        userfcn_fade(handles,'method1',0);
    case 'radiobutton_2'
        userfcn_fade(handles,'method2',0);
    otherwise
        disp('User-defined Warning: No radio button tag is ',...
           'found when SelectionChangeFcn is triggered.');
        return;
end

%updates the handles structure
guidata(hObject, handles);

function userfcn_fade(handles, method, onoff)
% Hightlight or fade selected method contents on input panel
userData = get(handles.figure1,'UserData');
switch method
    case 'method1'
        if onoff
            set(handles.listbox, 'Enable','on','BackgroundColor','white');
            set(handles.pushbutton_add,'Enable','on');
            set(handles.pushbutton_del,'Enable','on');
            set(handles.uipanel_2,'ForegroundColor','black');
            set(handles.edit_ps,'Enable','on');
            set(handles.edit_ti,'Enable','on');
            set(handles.text_body2,'ForegroundColor','black');
            set(handles.text_body3,'ForegroundColor','black');
            set(handles.pushbutton_done,'string','Save & Apply');
            % Save and reset text in Description text box
            setappdata(handles.edit_notes,'backup2',...
                get(handles.edit_notes,'String'));            
            set(handles.edit_notes,'String',...
                getappdata(handles.edit_notes,'backup1'));
            set(handles.pushbutton_output, 'Enable', 'on')
        else
            set(handles.listbox, 'Enable','off', 'BackgroundColor', [.95 .95 .95]);
            set(handles.pushbutton_add,'Enable','off');
            set(handles.pushbutton_del,'Enable','off');
            set(handles.uipanel_2,'ForegroundColor',[.5 .5 .5]);
            set(handles.edit_ps,'Enable','off');
            set(handles.edit_ti,'Enable','off');
            set(handles.text_body2,'ForegroundColor',[.5 .5 .5]);
            set(handles.text_body3,'ForegroundColor',[.5 .5 .5]);
            set(handles.pushbutton_output, 'Enable', 'off')
        end
    case 'method2'
        if onoff
            set(handles.pushbutton_open,'Enable','on');
            set(handles.pushbutton_done,'string','Apply');
            % Save and reset text in Description text box
            setappdata(handles.edit_notes,'backup1',...
                get(handles.edit_notes,'String'));
            set(handles.edit_notes,'String',...
                getappdata(handles.edit_notes,'backup2'));
            if isempty(userData.MD)
                set(handles.text_body8,'Visible','off');
            else
                set(handles.text_body8,'Visible','on');
            end
            
        else
            set(handles.pushbutton_open,'Enable','off'); 
            
            set(handles.text_body8,'Visible','off');

        end
    otherwise
        return;
end


% --------------------------------------------------------------------
function menu_file_Callback(hObject, eventdata, handles)
% hObject    handle to menu_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_about_Callback(hObject, eventdata, handles)
% hObject    handle to menu_about (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_about_lccb_Callback(hObject, eventdata, handles)
% hObject    handle to menu_about_lccb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
status = web('http://lccb.hms.harvard.edu/', '-browser');
if status
    switch status
        case 1
            msg = 'System default web browser is not found.';
        case 2
            msg = 'System default web browser is found but could not be launched.';
        otherwise
            msg = 'Fail to open browser for unknown reason.';
    end
    warndlg(msg,'Fail to open browser','modal');
end



% --------------------------------------------------------------------
function menu_file_open_Callback(hObject, eventdata, handles)
% Call back function of Manu File-Open

% Prepare eventdata for SelectionChangeFcn of radio panel
eventdata2 = struct('EventName','SelectionChanged',...
    'OldValue',get(handles.uipanel_1, 'SelectedObject'),...
    'NewValue', handles.radiobutton_2);
% Set radiobutton method2 selected
set(handles.uipanel_1, 'SelectedObject', handles.radiobutton_2);
% Go to SelectionChangeFcn for 
uipanel_1_SelectionChangeFcn(handles.uipanel_1, eventdata2);

%
pushbutton_open_Callback(handles.pushbutton_open,[ ],handles);



% --------------------------------------------------------------------
function menu_file_quit_Callback(hObject, eventdata, handles)
% hObject    handle to menu_file_quit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
delete(handles.figure1);




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


% --- Executes on selection change in listbox_2.
function listbox_2_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns listbox_2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_2


% --- Executes during object creation, after setting all properties.
function listbox_2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_ps2_Callback(hObject, eventdata, handles)
% hObject    handle to edit_ps2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_ps2 as text
%        str2double(get(hObject,'String')) returns contents of edit_ps2 as a double


% --- Executes during object creation, after setting all properties.
function edit_ps2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_ps2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_ti2_Callback(hObject, eventdata, handles)
% hObject    handle to edit_ti2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_ti2 as text
%        str2double(get(hObject,'String')) returns contents of edit_ti2 as a double


% --- Executes during object creation, after setting all properties.
function edit_ti2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_ti2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function menu_file_new_Callback(hObject, eventdata, handles)
% hObject    handle to menu_file_new (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Call back function of Menu File-New

% Prepare eventdata for SelectionChangeFcn of radio panel
eventdata2 = struct('EventName','SelectionChanged',...
    'OldValue',get(handles.uipanel_1, 'SelectedObject'),...
    'NewValue', handles.radiobutton_1);
% Set radiobutton method2 selected
set(handles.uipanel_1, 'SelectedObject', handles.radiobutton_1);
% Go to SelectionChangeFcn for 
uipanel_1_SelectionChangeFcn(handles.uipanel_1, eventdata2);

% --- Executes on button press in pushbutton_cancel.
function pushbutton_cancel_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
delete(handles.figure1);


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
delete(hObject);


% --- Executes during object deletion, before destroying properties.
function figure1_DeleteFcn(hObject, eventdata, handles)
userData = get(handles.figure1, 'UserData');

if isfield(userData, 'mainFig') && ishandle(userData.mainFig)
   setappdata(userData.mainFig, 'setupMovieDataFlag', 0); 
end

if isfield(userData, 'iconHelpFig') && ishandle(userData.iconHelpFig)
   delete(userData.iconHelpFig) 
end


% --- Executes on button press in checkbox_continue.
function checkbox_continue_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_continue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_continue



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


% --------------------------------------------------------------------
function menu_about_update_Callback(hObject, eventdata, handles)
status = web('http://lccb.hms.harvard.edu/software.html', '-browser');
if status
    switch status
        case 1
            msg = 'System default web browser is not found.';
        case 2
            msg = 'System default web browser is found but could not be launched.';
        otherwise
            msg = 'Fail to open browser for unknown reason.';
    end
    warndlg(msg,'Fail to open browser','modal');
end
