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

% Last Modified by GUIDE v2.5 23-Apr-2010 09:42:07

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
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to setupMovieDataGUI (see VARARGIN)

% Choose default command line output for setupMovieDataGUI
handles.output = hObject;

% Set callback function of radio button group uipanel_1
 set(handles.uipanel_1, 'SelectionChangeFcn', @uipanel_1_SelectionChangeFcn);
% Set radio button unchecked 
 set(handles.uipanel_1, 'SelectedObject', []);
 
% Set GUI data used to save MovieData object
handles.MD = [ ];

% Indicate the first package control panel the GUI will go to after 
% MovieData is sucessfully created and MovieDate setup panel is destroied.
% For future needs of multiple packages or multiple movie data input. 
% These variables can be re-defined as input variables when MovieData 
% setup panel is created.
handles.firstPackageName = 'BioSensorsPackage';
handles.firstPackageCtr = @BioSensorsPackage;

% Load help icon from dialogicons.mat
load lccbGuiIcons.mat
supermap(1,:) = get(hObject,'color');

axes(handles.axes_1);
Img = image(questIconData);
set(hObject,'colormap',supermap);
set(gca, 'XLim',get(Img,'XData'),'YLim',get(Img,'YData'),...
    'visible','off');
set(Img,'ButtonDownFcn',@help_ButtonDownFcn);

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

if isempty (get(handles.uipanel_1, 'SelectedObject'));
    errordlg('Please select a method to choose input channels',...
        'Please select a method','modal');
    return;
end
    
switch get(get(handles.uipanel_1, 'SelectedObject'), 'tag')
    case 'radiobutton_1'
        % Get input parameters channel path (cell array), pixel size and 
        % time interval
        channelPath = get(handles.listbox, 'string');
        pixelSize = get(handles.edit_ps, 'string');
        timeInterval = get(handles.edit_ti, 'string');
        notes = get(handles.edit_notes, 'string');
        % Check if text box is empty
        if isempty(channelPath)
            % Exception
            errordlg('Please provide at least one channel path',...
                       'User Input Error','modal');
            return;
        end
        if isempty(pixelSize) || ...
                isempty(timeInterval)
            % Exception
            errordlg('Please provide pixel size and time interval values',...
                'User Input Error','modal');
            return;
        end
        % Ask the user where to save the movie data
        [file,path] = uiputfile('*.mat','Save Movie Infomation as',...
            'movieData.mat');
        if ~any([file,path])
            return;
        end
        % Creat the movieData object (try/catch block) -------------!!!
        try
        MD = MovieData(channelPath,pixelSize,timeInterval,path,file,notes);
        catch ME
            errordlg([ME.message ' Input data is not saved.'],...
                'User Input Error','modal');
            return;
        end   

        try
        % Sanity Check
        if MD.sanityCheck(path, file, false);
            save([path file], 'MD');
            disp('Input data has been saved');
            % Save MovieDate as GUI data
             handles.MD = MD;
            guidata(hObject, handles);
        else 
            disp('User-defined: Potential problem exists in MovieData object');
        end
        
        % Catch: exception occurs - data is not saved.
        catch ME
            delete(MD);
            handles.MD = [ ];
            errordlg([ME.message 'Input data is not saved.'],...
                'Channel Path & Movie Data Error','modal');
            return;
        end

    case 'radiobutton_2'
        % Sanity check will be done in callback fcn of pushbutton_open 
        % Go to the first package
        if isempty(handles.MD)
            errordlg(['Movie data has not been sucessfully loaded. Please' ...
                ' click ''Open File...'' button to select a valid MAT file.'],...
                'Movie Data Not Loaded','modal');
            return;
        end
        % handles to the same MovieData object
        MD = handles.MD;
    otherwise
        error('User-defined: unexpected radio button has been selected');
end

% Go to the first package, call the control panel window
% handle.MD and MD are handles of the same MovieData object. Use MD for
% convenience
disp('Continue to the first package');
% I. Befor going to the first package, firstly check if the first package
% already exists    
packageExist = false;
if isempty(MD.packages_)
    % create a new package object
    MD.addPackage( handles.firstPackageCtr(MD) );
    % Set the current package to the newly created package
    MD.setCrtPackage( MD.packages_{1} ) ;
else
    % find existing package

    for i = 1: length(MD.packages_)
        if strcmp(class(MD.packages_{i}), handles.firstPackageName)
            % TODO: dlg box: previous package exists; Save package and
            % proceeds  
            
            MD.setCrtPackage(MD.packages_{i});
            packageExist = true;
            break;
        end
    end
    if ~packageExist
        % No same package is found.. create a new package object
        MD.addPackage( handles.firstPackageCtr(MD) )
        MD.setCrtPackage( MD.packages_{end} );
    end
end

% II. Before going to the first package, secondly check if existing
% processes can be recycled. New package has no processes in it.
if ~packageExist && ~isempty(MD.processes_)
    for i = 1: length(MD.crtPackage_.processClassNames_)
        % ith process in package's process list
        % 'sameProcID' save the ID of existing processes in MovieData that 
        % could be recycled by current package. 

        sameProcID = [];
        sameProcName = {};
        for j = 1: length(MD.processes_)
            % jth available process in MovieData's process pool
           if strcmp(class(MD.processes_{j}), ...
                        MD.crtPackage_.processClassNames_{i}) 
                sameProcID = horzcat (sameProcID, j);
                sameProcName = horzcat(sameProcName, ...
                                            {MD.processes_{j}.name_});
           end
        end
        if ~isempty(sameProcID)
            [select, ok] = listdlg('ListString', sameProcName, ...
                 'SelectionMode','single','ListSize',[420 200], ...
                 'Name','Select an existing process','PromptString', ...
                 {['There are existing ',sameProcName{1},' processes saved during previous sessions. You can choose to:'],...
                 '',...
                 '1. Use one of the following existing processes by click ''Select''',...
                 '2. click ''New'' to start over this process'} , ...
                 'OKString','Select','CancelString','New'); 
             if ok 
                 % Add process to current package
                 MD.crtPackage_.setProcess(i,MD.processes_{select});
             end
        end
    end
end


guidata(hObject, handles);


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
% hObject    handle to pushbutton_add (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
path = uigetdir;
if path == 0
    return;
end
% Input validation function ...
% Get current list
contents = get(handles.listbox,'String');
if any(strcmp(contents,path))
   warndlg('Entry channel exists! Please select a differenct directory.',...
       'Warning','modal');
   return; 
end
% Add current formula to the listbox
contents{end+1} = path;
set(handles.listbox,'string',contents);

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
% hObject    handle to pushbutton_open (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uigetfile('*.mat','Select one MAT-file');
if ~any([filename,pathname])
    return;
end
try
    pre = whos('-file', [pathname filename]);  % - Exception: fail to access .mat file
catch ME
    errordlg(ME.message,'Fail to open MAT file','modal');
    return;
end
% Find MovieData object in .mat file before loading it
structMD = pre( logical(strcmp({pre(:).class},'MovieData')) );
switch length(structMD)
    case 0
        % Exception: No MovieData object is in the MAT file
        errordlg('No MovieData object is found in selected .mat file',...
            'Wrong .MAT File','modal');
        return;
    case 1
        % Load MovieData object in .mat file
        load([pathname filename],'-mat',structMD.name);
        % Name/Rename MovieData obj as 'MD'
        eval(['MD =' structMD.name ';']);
    otherwise
        % Exception - TODO: multiple movies saved in one MAT file
        errordlg('More than one movies are saved in .mat file',...
            'Bad .MAT File','modal');
        return;
end
% Sanity check of saved movie data
try
    % Set parameter values in .mat file to GUI 
        set(handles.edit_mat,'String',[pathname,filename]);
        set(handles.listbox_2,'String',MD.channelPath_);
        set(handles.edit_ps2,'String',MD.pixelSize_);
        set(handles.edit_ti2,'String',MD.timeInterval_);
        set(handles.edit_notes,'String',MD.notes_)
    
    if MD.sanityCheck( pathname,filename,false );       
        % Success Notice: Movie Data has been loaded
        set(handles.text_body8,'visible','on');
        handles.MD = MD;
        guidata(hObject,handles);
    end
catch ME
    % If don't pass sanity check, set MovieData unloaded
    delete(MD);
    handles.MD = [ ];
    errordlg([ME.message 'Movie data is not sucessfully loaded.'],...
                'Channel Path & Movie Data Error','modal');   
    % Erase success notice
    set(handles.text_body8,'Visible','off');
    guidata(hObject, handles);
    
    return;
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
    % By the way, Set notes text box enabled when a radio button is 
    % selected for the first time
    
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
            set(handles.pushbutton_done,'string','Save & Finish');
            set(handles.edit_notes, 'enable', 'on');
            % Save and reset text in Description text box
            setappdata(handles.edit_notes,'backup2',...
                get(handles.edit_notes,'String'));            
            set(handles.edit_notes,'String',...
                getappdata(handles.edit_notes,'backup1'));
        else
            set(handles.listbox, 'Enable','off', 'BackgroundColor', [.95 .95 .95]);
            set(handles.pushbutton_add,'Enable','off');
            set(handles.pushbutton_del,'Enable','off');
            set(handles.uipanel_2,'ForegroundColor',[.5 .5 .5]);
            set(handles.edit_ps,'Enable','off');
            set(handles.edit_ti,'Enable','off');
            set(handles.text_body2,'ForegroundColor',[.5 .5 .5]);
            set(handles.text_body3,'ForegroundColor',[.5 .5 .5]);
            set(handles.edit_notes, 'enable', 'off');

        end
    case 'method2'
        if onoff
            set(handles.pushbutton_open,'Enable','on');
            set(handles.pushbutton_done,'string','Finish');
            % Save and reset text in Description text box
            setappdata(handles.edit_notes,'backup1',...
                get(handles.edit_notes,'String'));
            set(handles.edit_notes,'String',...
                getappdata(handles.edit_notes,'backup2'));
            if isempty(handles.MD)
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
function menu_help_lccb_Callback(hObject, eventdata, handles)
% hObject    handle to menu_help_lccb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


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

function help_ButtonDownFcn(hObject, eventdata)
% Call back function when help icon is clicked
t = help('plot');
helpdlg(t,'Help');
