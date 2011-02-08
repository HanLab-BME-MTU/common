function varargout = uTrackPackageGUI(varargin)
% UTRACKPACKAGEGUI M-file for uTrackPackageGUI.fig
%      UTRACKPACKAGEGUI, by itself, creates a new UTRACKPACKAGEGUI or raises the existing
%      singleton*.
%
%      H = UTRACKPACKAGEGUI returns the handle to a new UTRACKPACKAGEGUI or the handle to
%      the existing singleton*.
%
%      UTRACKPACKAGEGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in UTRACKPACKAGEGUI.M with the given input arguments.
%
%      UTRACKPACKAGEGUI('Property','Value',...) creates a new UTRACKPACKAGEGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before uTrackPackageGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to uTrackPackageGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help uTrackPackageGUI

% Last Modified by GUIDE v2.5 28-Nov-2010 15:51:53

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @uTrackPackageGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @uTrackPackageGUI_OutputFcn, ...
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


% --- Executes just before uTrackPackageGUI is made visible.
function uTrackPackageGUI_OpeningFcn(hObject, eventdata, handles, varargin)
%
% uTrackPackageGUI(MD)   MD: MovieData object
%
% Useful tools
%
% User Data:
%
%       userData.MD - array of MovieData object
%       userData.package - array of package (same length with userData.MD)
%       userData.crtPackage - the package of current MD
%       userData.id - the id of current MD on board
%
%       userData.dependM - dependency matrix
%       userdata.statusM - GUI status matrix
%       userData.optProcID - optional process ID
%
%       userData.passIconData - pass icon image data
%       userData.errorIconData - error icon image data
%       userData.warnIconData - warning icon image data
%       userData.questIconData - help icon image data
%       userData.colormap - color map
%
%       userData.setFig - array of handles of (multiple) setting figures (may not exist)
%       userData.resultFig - array of handles of (multiple) result figures (may not exist)
%       userData.setupMovieDataFig - handle of (single) setupMovieData figure (may not exist)
%       userData.overviewFig - handles of (single) overviewMovieDataGUI figure (may not exist)
%       userData.packageHelpFig - handle of (single) help figure (may not exist)
%       userData.iconHelpFig - handle of (single) help figures (may not exist)
%       userData.processHelpFig - handle of (multiple) help figures (may not exist) 
%       userData.msgboxGUI - handle of message box
%       
%
% NOTE:
%   
%   userData.statusM - 1 x m stucture array, m is the number of Movie Data 
%                      this user data is used to save the status of movies
%                      when GUI is switching between different movie(s)
%                   
%   	fields: IconType - the type of status icons, 'pass', 'warn', 'error'
%               Msg - the message displayed when clicking status icons
%               Checked - 1 x n logical array, n is the number of processes
%                         used to save value of check box of each process
%               Visited - logical true or false, if the movie has been
%                         loaded to GUI before 
%


packageName = 'UTrackPackage';

% Load movie data and recycle processes
userfcn_iniPackage_commonCode;


% --- Outputs from this function are returned to the command line.
function varargout = uTrackPackageGUI_OutputFcn(hObject, eventdata, handles) 
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
userData = get(handles.figure1, 'UserData');
for i = 1: length(userData.MD)
    userData.MD(i).saveMovieData
end
delete(handles.figure1);

% --- Executes on button press in pushbutton_status.
function pushbutton_status_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_status (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

userData = get(handles.figure1, 'UserData');

% if newMovieDataGUI exist
if isfield(userData, 'overviewFig') && ishandle(userData.overviewFig)
    delete(userData.overviewFig)
end

userData.overviewFig = newMovieDataGUI('mainFig',handles.figure1, 'overview', userData.MD(userData.id));
set(handles.figure1, 'UserData', userData);

% --- Executes on button press in pushbutton_save.
function pushbutton_save_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

userData = get(handles.figure1, 'UserData');

for i = 1: length(userData.MD)
    userData.MD(i).saveMovieData
end

set(handles.text_body3, 'Visible', 'on')
pause(1)
set(handles.text_body3, 'Visible', 'off')

% --- Executes on button press in pushbutton_left.
function pushbutton_left_Callback(hObject, eventdata, handles)
% userData.id
% userData.crtPackage
%
userData = get(handles.figure1, 'UserData');
l = length(userData.MD);

userData.statusM(userData.id).Checked = userfcn_saveCheckbox(handles);

userData.id = userData.id - 1;

if userData.id < 1
   userData.id = l;
end

userData.crtPackage = userData.package(userData.id);
set(handles.figure1, 'UserData', userData)

% Set up movie explorer
set(handles.popupmenu_movie, 'Value', userData.id)

% Set up GUI
if userData.statusM(userData.id).Visited
   userfcn_updateGUI(handles, 'refresh') 
else
   userfcn_updateGUI(handles, 'initialize') 
end

% --- Executes on button press in pushbutton_right.
function pushbutton_right_Callback(hObject, eventdata, handles)
userData = get(handles.figure1, 'UserData');
l = length(userData.MD);

userData.statusM(userData.id).Checked = userfcn_saveCheckbox(handles);

userData.id = userData.id + 1;

if userData.id > l
   userData.id = mod(userData.id, l);
end

userData.crtPackage = userData.package(userData.id);
set(handles.figure1, 'UserData', userData)

% Set up movie explorer
set(handles.popupmenu_movie, 'Value', userData.id)

% Set up GUI
if userData.statusM(userData.id).Visited
   userfcn_updateGUI(handles, 'refresh') 
else
   userfcn_updateGUI(handles, 'initialize') 
end


% --- Executes on selection change in popupmenu_movie.
function popupmenu_movie_Callback(hObject, eventdata, handles)

userData = get(handles.figure1, 'UserData');

if get(hObject, 'Value') == userData.id
   return 
end

l = length(userData.MD);
userData.statusM(userData.id).Checked = userfcn_saveCheckbox(handles);

userData.id = get(hObject, 'Value');
userData.crtPackage = userData.package(userData.id);
set(handles.figure1, 'UserData', userData)

% Set up GUI
if userData.statusM(userData.id).Visited
   userfcn_updateGUI(handles, 'refresh') 
else
   userfcn_updateGUI(handles, 'initialize') 
end



% --- Executes during object creation, after setting all properties.
function popupmenu_movie_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_movie (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
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


% --- Executes on button press in checkbox_1.
function checkbox_1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_1
userfcn_checkAllMovies(1, get(hObject,'value'), handles);
userfcn_lampSwitch(1, get(hObject,'value'), handles);


% --- Executes on button press in pushbutton_set_1.
function pushbutton_set_1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_set_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
userData = get(handles.figure1, 'UserData');
procID = 1;
userData.setFig(procID) = detectionProcessGUI('mainFig',handles.figure1,procID);
set(handles.figure1, 'UserData', userData);
guidata(hObject,handles);

% --- Executes on button press in pushbutton_show_1.
function pushbutton_show_1_Callback(hObject, eventdata, handles)

procID = 1;
userData = get(handles.figure1, 'UserData');

% ----------------------------------
crtProc = userData.crtPackage.processes_{procID};
% Make sure output exists
chan = [];
for i = 1:length(userData.MD(userData.id).channels_)
    
    if ~isempty(crtProc.outParams_{i})
        chan = i; 
        break
    end
end

if isempty(chan)
   warndlg('The current step does not have any output yet.','No Output','modal');
   return
end

% Make sure detection output is valid
firstframe = [];
for i = 1:length(crtProc.outParams_{chan}.movieInfo)
   
    if ~isempty(crtProc.outParams_{chan}.movieInfo(i).amp)
        firstframe = i;
        break
    end
end

if isempty(firstframe)
    warndlg('The detection result is empty. There is nothing to visualize.','Empty Output','modal');
   return
end
% -------------------------------------

if isfield(userData, 'resultFig') && ishandle(userData.resultFig)
    
    delete(userData.resultDet)
end
%     userData.resultDet = userData.crtPackage.processes_{procID}.showResult;
    userData.resultDet = detectionVisualGUI('mainFig', handles.figure1, procID);


set(handles.figure1, 'UserData', userData);


% --- Executes on button press in checkbox_2.
function checkbox_2_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_2
userfcn_checkAllMovies(2, get(hObject,'value'), handles);
userfcn_lampSwitch(2, get(hObject,'value'), handles);

% --- Executes on button press in pushbutton_set_2.
function pushbutton_set_2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_set_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
userData = get(handles.figure1, 'UserData');
procID = 2;
% userData.setFig(procID) = segmentProcessGUI('mainFig',handles.figure1,procID);
userData.setFig(procID) = trackingProcessGUI('mainFig',handles.figure1,procID);
set(handles.figure1, 'UserData', userData);
guidata(hObject,handles);

% --- Executes on button press in pushbutton_show_2.
function pushbutton_show_2_Callback(hObject, eventdata, handles)

procID = 2;
userData = get(handles.figure1, 'UserData');

% ----------------------------------------------
crtProc = userData.crtPackage.processes_{procID};
% Make sure output exists
chan = [];
for i = 1:length(userData.MD(userData.id).channels_)
    
    if ~isempty(crtProc.outParams_{i})
        chan = i; 
        break
    end
end

if isempty(chan)
   warndlg('The current step does not have any output yet.','No Output','modal');
   return
end

% Make sure detection output is valid

if isempty(crtProc.outParams_{chan}.tracksFinal)
    warndlg('The tracking result is empty. There is nothing to visualize.','Empty Output','modal');
    return
end

% --------------------------------------------------

if isfield(userData, 'resultFig') && ishandle(userData.resultFig)
    
    delete(userData.resultDet)
end
%     userData.resultDet = userData.crtPackage.processes_{procID}.showResult;
    userData.resultDet = trackingVisualGUI('mainFig', handles.figure1, procID);


set(handles.figure1, 'UserData', userData);


% --- Executes on button press in pushbutton_run.
function pushbutton_run_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% userfcn_pushbutton_run_common


% This is a common section of code called by pushbutton_run_Callback
% when user click the "Run" button on package control panels.
%
% Chuangang Ren
% 11/2010

userData = get(handles.figure1,'UserData');

movieRun = []; % id of movie to run
numMD = length(userData.MD); % number of movies


% Get check box status of current movie and update user data
userData.statusM(userData.id).Checked = userfcn_saveCheckbox(handles);
set(handles.figure1, 'UserData', userData)

% Determine the movie(s) to be processed
if ~get(handles.checkbox_runall, 'Value')
    
    if any(userData.statusM(userData.id).Checked)
       movieRun = cat(2, movieRun, userData.id);
    end
else
    
    % Check other movies
    for i = 0:numMD-1
        k = mod( userData.id + i, numMD);
        if k == 0 
            k = numMD; 
        end
    
        if any(userData.statusM(k).Checked)
            movieRun = cat(2, movieRun, k);
        end
    end     
end

if isempty(movieRun)
    warndlg('No step is selected, please select a step to process.','No Step Selected','modal');
    return
end



% ----------------------- Pre-processing examination ----------------------


% movie exception (same length of movie data)
movieException = cell(1, numMD);

procCheck = cell(1, numMD);%  id of checked processes 
procRun = cell(1, numMD);%  id of processes to run
optionalProcID = cell(1, numMD);% id of first-time-run optional process

for x = movieRun
    

procCheck{x} = find(userData.statusM(x).Checked); 

% Check if process exist
for i = procCheck{x}
    if isempty (userData.package(x).processes_{i})
        
        ME = MException('lccb:run:setup', 'Step %d is not set up yet. Tip: when step is set up successfully, the step name becomes bold.',i);
        movieException{x} = cat(2, movieException{x}, ME);
        
    end    
end

if ~isempty(  movieException{x} )
   continue 
end

% Check if selected processes have alrady be successfully run
% If force run, re-run every process that is checked
if ~get(handles.checkbox_forcerun, 'Value')

    k = true;
    for i = procCheck{x}

        if  ~( userData.package(x).processes_{i}.success_ && ...
            ~userData.package(x).processes_{i}.procChanged_ ) || ...
            ~userData.package(x).processes_{i}.updated_
        
            k = false;
            procRun{x} = cat(2, procRun{x}, i);
        end
    end
    if k
        movieRun = setdiff(movieRun, x);
        continue
    end
else
    procRun{x} = procCheck{x};
end



% Package full sanity check. Sanitycheck every checked process
procEx = userData.package(x).sanityCheck(true, procRun{x});

% Return user data !!!
set(handles.figure1, 'UserData', userData)

for i = procRun{x}
   if ~isempty(procEx{i})
       
       % Check if there is fatal error in exception array
       if strcmp(procEx{i}(1).identifier, 'lccb:set:fatal') || ...
               strcmp(procEx{i}(1).identifier, 'lccb:input:fatal')
           
           % Sanity check error - switch GUI to the x th movie 
           if x ~= userData.id
             set(handles.popupmenu_movie, 'Value', x)
             popupmenu_movie_Callback(handles.popupmenu_movie, [], handles) % user data retrieved, updated and submitted
           end
           
           userfcn_drawIcon(handles,'error', i, procEx{i}(1).message, true); % user data is retrieved, updated and submitted
           
           ME = MException('lccb:run:sanitycheck', 'Step %d %s: \n%s', i,userData.package(x).processes_{i}.name_, procEx{i}(1).message);
           movieException{x} = cat(2, movieException{x}, ME);
           
       end
   end
end

% Refresh user data !!!
userData = get(handles.figure1, 'UserData');


end

% --------------------- pre-processing examination ends -------------------

% Ok, now all evils are in movieException (1 x movielength  cell array), if there is any
% if yes - abort program and popup a error report
% if no - continue to process movie data
if isempty(movieRun)
    warndlg('All selected steps have been processed successfully. Please check the ''Force Run'' check box if you want to re-process the successful steps.','No Step Selected','modal');
    return
end

temp = find(~cellfun(@(x)isempty(x), movieException, 'UniformOutput', true));

if ~isempty(temp)
    msg = [];
    for i = 1:length(temp)
        if i == 1
            msg = strcat(msg, sprintf('Movie %d - %s:', temp(i), userData.MD(temp(i)).movieDataFileName_));
        else
            msg = strcat(msg, sprintf('\n\n\nMovie %d - %s:', temp(i), userData.MD(temp(i)).movieDataFileName_));
        end
        for j = 1:length(movieException{temp(i)})
            msg = strcat(msg, sprintf('\n-- %s', movieException{temp(i)}(j).message));
        end

    end
    msg = strcat(msg, sprintf('\n\n\nPlease solve the above problems before continuing. The Movie(s) couldn’t be processed.'));
    titlemsg = sprintf('Processing could not be continued for the following reasons:');
    
    % if msgboxGUI exist
    if isfield(userData, 'msgboxGUI') && ishandle(userData.msgboxGUI)
        delete(userData.msgboxGUI)
    end    
    
    userData.msgboxGUI = msgboxGUI('title',titlemsg,'text', msg);
    return
    
end

% ------------------------ Start Processing -------------------------------
kk = 0;
for x = movieRun
    
kk = kk+1;    
if x ~= userData.id
    
    set(handles.popupmenu_movie, 'Value', x)
    set(handles.figure1, 'UserData', userData)
    
    popupmenu_movie_Callback(handles.popupmenu_movie, [], handles) % user data retrieved, updated and submitted
    userData = get(handles.figure1, 'UserData');
end
    
% Find first-time-run optional process ID
for i = intersect(procRun{x}, userData.optProcID);
    if ~userData.package(x).processes_{i}.success_
        optionalProcID{x} = cat(2, optionalProcID{x}, i);
    end
end

% Set all running processes' sucess = false; 
for i = procRun{x}
    userData.crtPackage.processes_{i}.setSuccess(false);
end

% Clear icons of selected processes
% Return user data !!!
set(handles.figure1, 'UserData', userData)
userfcn_drawIcon(handles,'clear',procRun{x},'',true); % user data is retrieved, updated and submitted
% Refresh user data !!!
userData = get(handles.figure1, 'UserData');

% Disable 'Run' button
set(handles.pushbutton_run, 'Enable', 'off')
set(handles.checkbox_forcerun, 'Enable', 'off')
set(handles.checkbox_runall, 'Enable', 'off')
set(handles.text_status, 'Visible', 'on')

% Run algorithms!
try
    % Return user data !!!
    set(handles.figure1, 'UserData', userData)
    
    for i = procRun{x}
        set(handles.text_status, 'String', sprintf('Step %d - Processing %d of totally %d movies ...', i, kk, length(movieRun)) )
        userfcn_runProc_dfs(i, procRun{x}, handles); % user data is retrieved, updated and submitted

    end
    
catch ME
    
    set(handles.pushbutton_run, 'Enable', 'on') %%%%%
    set(handles.checkbox_forcerun, 'Enable', 'on') %%%%%
    set(handles.checkbox_runall, 'Enable', 'on') %%%%%
    set(handles.text_status, 'Visible', 'off') %%%%%
    throw(ME) %%%%%
    
    % Save the error into movie Exception cell array
    movieException{x} = ME;
    
    procRun{x} = procRun{x}(procRun{x} < i);
    optionalProcID{x} = optionalProcID{x}(optionalProcID{x} < i);
    
    
end

% Refresh user data !!!
userData = get(handles.figure1, 'UserData');
set(handles.pushbutton_run, 'Enable', 'on')
set(handles.checkbox_forcerun, 'Enable', 'on')
set(handles.checkbox_runall, 'Enable', 'on')
set(handles.text_status, 'Visible', 'off')



% ------- Check optional processes ----------

% Return user data !!!
set(handles.figure1, 'UserData', userData)
% In here, optionalProcID are successfuly first-time-run optional process ID
if ~isempty(optionalProcID{x})
    
    procEx = userData.crtPackage.checkOptionalProcess(procRun{x}, optionalProcID{x});
    
    for i = 1:size(userData.dependM, 1)
        if ~isempty(procEx{i})
            
            userfcn_drawIcon(handles,'warn',i,procEx{i}(1).message, true); % user data is retrieved, updated and submitted
        
        end
    end
end

end

% ----------------------------- Create error report ---------------------------------------


temp = find(~cellfun(@(x)isempty(x), movieException, 'UniformOutput', true));

if ~isempty(temp)
    msg = [];
    for i = 1:length(temp)
        if i == 1
            msg = strcat(msg, sprintf('Movie %d - %s:\n\n%s', ...
                temp(i), userData.MD(temp(i)).movieDataFileName_, movieException{x}.message));
        else
            msg = strcat(msg, sprintf('\n\n\nMovie %d - %s:\n\n%s', ...
                temp(i), userData.MD(temp(i)).movieDataFileName_, movieException{x}.message));
        end
    end
   
    msg = strcat(msg, sprintf('\n\n\nPlease verify your settings are correct. Feel free to contact us if you have question regarding this error.\n\nPlease help us improve the software by clearly reporting the scenario when this error occurs, and the above error information to us (error information is also displayed in Matlab command line).\nFor contact information please refer to the following URL:\n\nlccb.hms.harvard.edu/software.html'));

    % if msgboxGUI exist
    if isfield(userData, 'msgboxGUI') && ishandle(userData.msgboxGUI)
        delete(userData.msgboxGUI)
    end   
    if length(temp) == 1
        userData.msgboxGUI = msgboxGUI('title','The processing of following movie is terminated by run time error:','text', msg); 
    else
        userData.msgboxGUI = msgboxGUI('title','The processing of following movies are terminated by run time errors:','text', msg); 
    end
    
elseif length(movieRun) > 1 
    userData.iconHelpFig = helpdlg('All your movie data has been processed successfully.', 'Biosensors Package');
    set(handles.figure1, 'UserData', userData)
end


% --- Executes on button press in checkbox_forcerun.
function checkbox_forcerun_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_forcerun (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_forcerun


% --- Executes on button press in checkbox_runall.
function checkbox_runall_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_runall (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_runall


% --- Executes on button press in checkbox_all.
function checkbox_all_Callback(hObject, eventdata, handles)

% userData = get(handles.figure1, 'UserData');
% l = size(userData.dependM, 1);

% switch get(hObject,'value')
%     case 1
%         userfcn_enable(1:length(userData.dependM(:,1)),'on',handles,true);
        
        % track checked status
%         for x = 1: length(userData.MD)
%             userData.statusM(x).Checked = ones(1, l);
%         end
        
%     case 0
%         k = [];
%         for i = 1: size(userData.dependM, 1)
%             if ~isempty(userData.crtPackage.processes_{i}) && ...
%                 userData.crtPackage.processes_{i}.success_ 
%                 k = [k, i];
%             end
%             eval( ['set(handles.checkbox_',num2str(i),',''value'',0)']  );
%         end
%         tempDependM = userData.dependM;
%         tempDependM(:,k) = zeros(size(userData.dependM,1),length(k));
%         userfcn_enable(find(any(tempDependM,2)),'off',handles,true);
%         
%         % track checked status
%         for x = 1: length(userData.MD)
%             userData.statusM(x).Checked = zeros(1, l);
%         end        
        
% end

% set(handles.figure1, 'UserData', userData)


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
userData = get(handles.figure1,'Userdata');
if isfield(userData, 'MD')
    MD = userData.MD;
else
    delete(handles.figure1);
    return;
end

user_response = questdlg('Do you want to save the current progress?', ...
    'BioSensors Package Control Panel');
switch lower(user_response)
    case 'yes'
        for i = 1: length(userData.MD)
            userData.MD(i).saveMovieData
        end
        delete(handles.figure1);
    case 'no'
        delete(handles.figure1);
    case 'cancel'
end


% --- Executes during object deletion, before destroying properties.
function figure1_DeleteFcn(hObject, eventdata, handles)

userData = get(handles.figure1, 'UserData');
% setFlag = getappdata(hObject, 'setFlag');

% Delete setting figures (multiple)
if isfield(userData, 'setFig')
    for i = 1: length(userData.setFig)
        if userData.setFig(i)~=0 && ishandle(userData.setFig(i))
            delete(userData.setFig(i))
        end
    end
end

% Close result figures (multiple)
if isfield(userData, 'resultFig') && userData.resultFig~=0 && ishandle(userData.resultFig)

	delete(userData.resultFig(i))
end

% If open, delete MovieData Overview GUI figure (single)
if isfield(userData, 'setupMovieDataFig') && ishandle(userData.setupMovieDataFig)
   delete(userData.setupMovieDataFig) 
end

% If open, delete MovieData Overview GUI figure (single)
if isfield(userData, 'overviewFig') && ishandle(userData.overviewFig)
   delete(userData.overviewFig) 
end

% Delete pre-defined package help dialog (single)
if isfield(userData, 'packageHelpFig') && ishandle(userData.packageHelpFig)
   delete(userData.packageHelpFig) 
end

% Delete pre-defined icon help dialog (single)
if isfield(userData, 'iconHelpFig') && ishandle(userData.iconHelpFig)
   delete(userData.iconHelpFig) 
end

% Delete pre-defined process help dialogssetting figures (multiple)
if isfield(userData, 'processHelpFig')
    for i = 1: length(userData.processHelpFig)
        if userData.processHelpFig(i)~=0 && ishandle(userData.processHelpFig(i))
            delete(userData.processHelpFig(i))
        end
    end
end

% msgboxGUI used for error reports
if isfield(userData, 'msgboxGUI') && ishandle(userData.msgboxGUI)
   delete(userData.msgboxGUI) 
end


% --- Executes on key press with focus on figure1 and none of its controls.
function figure1_KeyPressFcn(hObject, eventdata, handles)

if strcmp(eventdata.Key, 'return')
    pushbutton_done_Callback(handles.pushbutton_done, [], handles);
end
if strcmp(eventdata.Key, 'leftarrow')
    pushbutton_left_Callback(handles.pushbutton_left, [], handles);
end
if strcmp(eventdata.Key, 'rightarrow')
    pushbutton_right_Callback(handles.pushbutton_right, [], handles);
end


% --------------------------------------------------------------------
function menu_file_Callback(hObject, eventdata, handles)
% hObject    handle to menu_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_help_lccb_Callback(hObject, eventdata, handles)
% hObject    handle to menu_help_lccb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


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
% Call back function of 'New' in menu bar
userData = get(handles.figure1,'Userdata');

        for i = 1: length(userData.MD)
            userData.MD(i).saveMovieData
        end

movieSelectorGUI
delete(handles.figure1)


% --------------------------------------------------------------------
function menu_file_save_Callback(hObject, eventdata, handles)

userData = get(handles.figure1, 'UserData');
userData.MD(userData.id).saveMovieData

set(handles.text_body3, 'Visible', 'on')
pause(1)
set(handles.text_body3, 'Visible', 'off')


% --------------------------------------------------------------------
function menu_file_exit_Callback(hObject, eventdata, handles)
% hObject    handle to menu_file_exit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
delete(handles.figure1);
