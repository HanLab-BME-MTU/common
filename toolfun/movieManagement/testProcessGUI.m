function varargout = testProcessGUI(varargin)
% TESTPROCESSGUI M-file for testProcessGUI.fig
%      TESTPROCESSGUI, by itself, creates a new TESTPROCESSGUI or raises the existing
%      singleton*.
%
%      H = TESTPROCESSGUI returns the handle to a new TESTPROCESSGUI or the handle to
%      the existing singleton*.
%
%      TESTPROCESSGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TESTPROCESSGUI.M with the given input arguments.
%
%      TESTPROCESSGUI('Property','Value',...) creates a new TESTPROCESSGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before testProcessGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to testProcessGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help testProcessGUI

% Last Modified by GUIDE v2.5 04-May-2010 10:30:38

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @testProcessGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @testProcessGUI_OutputFcn, ...
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


% --- Executes just before testProcessGUI is made visible.
function testProcessGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% Available tools 
% GUI data:
%       handles.mainFig - handle of main figure
%       handles.handles_main - 'handles' of main figure
%       handles.procID - The ID of process in the current package
%       handles.crtProc - handle of current process
%       handles.procConstr - constructor of current process

userData = get(handles.figure1, 'UserData');

% Choose default command line output for testProcessGUI
handles.output = hObject;

t = find(strcmp(varargin,'mainFig'));
userData.mainFig = varargin{t+1};
userData.procID = varargin{t+2};
userData.handles_main = guidata(userData.mainFig);

% Get current package
userData_main = get(userData.mainFig, 'UserData');
userData.crtPackage = userData_main.crtPackage;

userData.crtProc = userData.crtPackage.processes_{userData.procID};

eval ( [ 'userData.procConstr = @', ...
    userData.crtPackage.processClassNames_{userData.procID},';']); % No need

% Write introduction text
set(handles.text_body1,'string',['Setting for process ', ...
   userData.crtPackage.processClassNames_{userData.procID} ]); % No need

% Set existing parameters, if there is any

if ~isempty(userData.crtProc)
    set(handles.edit_para, 'string', num2str(userData.crtProc.funParams_));
end

% Set the flag of setting panel as 1
setFlag = getappdata(userData.mainFig,'setFlag');
setFlag(userData.procID) = 1;
setappdata(userData.mainFig,'setFlag',setFlag);

% Update user data
set(userData.mainFig, 'UserData', userData_main);
set(handles.figure1, 'UserData', userData);
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = testProcessGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = hObject;



function edit_para_Callback(hObject, eventdata, handles)
% If process exists (not a new process), and if text box is changed, set
% procChanged_ to 'true'
userData = get(handles.figure1, 'UserData');

if ~isempty( userData.crtProc )
    userData.crtProc.setProcChanged(true);
end

set(handles.figure1, 'UserData', userData);
guidata(hObject, handles)

% --- Executes during object creation, after setting all properties.
function edit_para_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_para (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% function pushbutton_done_Callback(hObject, eventdata, handles)
%
% % Get the function name from the GUI component (listbox for instance)
% funName = get(hListBox, 'String');
%
% % Get the pathToResult from a browser
% pathToResult = uigetdir(...);
%
% if ~ischar(pathToResult) % the user presses cancel
%    errodlg(pathToResult needs to be reselected !);
% en
%
% % Get the alpha parameter
% alpha = get(hAlpha, 'Value');
%
% if alpha < 0
%   errordlg(alpha value is wrong !!!);
% end
%
% Up to here, every parameter needed by segmentationProcess are right.
%
% % Find in the movieData.packages(currPackage).listOfProcesses whether
% % segmentation process is here
% ...
% if exist(segmentationProcess)
%    segmentationProcess.setFuncName(funcName);
%    segmentationProcess.setPathToResult(pathToResult);
%    segmentationProcess.setAlpha(alpha);
%    segmentationProcess.hasChanged = true;
% else
%    % Create a new segmentationProcess
%    newProc = procConstr(MD,funcName, pathToResult, alpha);
%
%    % Append the new process to the current package process list.
%    movieData.packages(currPackage).listOfProcess  =
%    [movieData.packages(currPackage).listOfProcess, newProc];
% end
%
% END OF CALLBACK

% segmentationProcess.sanityCheck(boolean full)
%   check that any other field (example: file result, etc.) exist and
%   contains the expected data.
%
%   the parameter 'funName', 'pathToResult' and 'alpha' don't need to be
%   checked here because we know that there are valid (validation done in
%   the 'finish' callback above.

% --- Executes on button press in pushbutton_done.
function pushbutton_done_Callback(hObject, eventdata, handles)

userData = get(handles.figure1, 'UserData');
userData_main = get(userData.mainFig, 'UserData');
MD = userData_main.MD;

para = str2double(get(handles.edit_para, 'string'));
% Check if input is blank
if isnan(para)
   errordlg('Pamameter cannot be blank',...
                       'User Input Error','modal');
   return;
end

pObj = userData.procConstr(MD, 'maskPath', 'funName', para);

% Do process sanity check - only check pamameter error
try
    pObj.sanityCheck;
catch ME
    delete(pObj);
    errordlg([ME.message 'Please double check your settings'],...
                'Setting Error','modal');
    return;
end

% If process sanity check passes, set new parameters
if isempty( userData.crtProc )
    % Add new process to both process lists of MovieData and current package
    MD.addProcess(pObj);
    userData.crtPackage.setProcess(userData.procID,pObj);
    % Save MovieData
    save([MD.movieDataPath_ MD.movieDataFileName_], 'MD');
else
    % Change existing parameters to the current setting
   userData.crtProc.setPara(para) 
   save([MD.movieDataPath_ MD.movieDataFileName_], 'MD');
end

% Do sanity check - only check changed parameters
procEx = userData.crtPackage.sanityCheck(false,'all');

% Draw some bugs on the wall :)
for i = 1: length(procEx)
   if ~isempty(procEx{i})
       
       % Set the processChanged to true
%        userData.crtPackage.processes_{i}.setProcChanged(true);
       % Draw warning label on the i th process
      userfcn_drawIcon(userData.handles_main,'warn',i,procEx{i}(1).message);
      
       for j = 1: length(procEx{i})
           if strcmp(procEx{i}(j).identifier, 'lccb:set:fatal');
               error(['User-defined: lccb:set:fatal exception should not', ...
                   ' be here. Check userData.crtPackage.sanityCheck']);
           end
       end
   end
end

% Save user data
set(userData.mainFig, 'UserData', userData_main)
set(handles.figure1, 'UserData', userData);
guidata(hObject,handles);
delete(handles.figure1);



% --- Executes on button press in pushbutton_cancel.
function pushbutton_cancel_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
delete(handles.figure1);


% --- Executes during object deletion, before destroying properties.
function figure1_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
userData = get(handles.figure1, 'UserData');

setFlag = getappdata(userData.mainFig,'setFlag');
setFlag(userData.procID) = 0;
setappdata(userData.mainFig,'setFlag',setFlag);

set(handles.figure1, 'UserData', userData);
guidata(hObject,handles);




