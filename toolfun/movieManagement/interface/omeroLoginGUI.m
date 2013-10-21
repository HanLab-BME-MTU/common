function varargout = omeroLoginGUI(varargin)
% OMEROLOGINGUI MATLAB code for omeroLoginGUI.fig
%      OMEROLOGINGUI, by itself, creates a new OMEROLOGINGUI or raises the existing
%      singleton*.
%
%      H = OMEROLOGINGUI returns the handle to a new OMEROLOGINGUI or the handle to
%      the existing singleton*.
%
%      OMEROLOGINGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in OMEROLOGINGUI.M with the given input arguments.
%
%      OMEROLOGINGUI('Property','Value',...) creates a new OMEROLOGINGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before omeroLoginGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to omeroLoginGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help omeroLoginGUI

% Last Modified by GUIDE v2.5 20-Oct-2013 22:06:57

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @omeroLoginGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @omeroLoginGUI_OutputFcn, ...
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


% --- Executes just before omeroLoginGUI is made visible.
function omeroLoginGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to omeroLoginGUI (see VARARGIN)


global client
global session

if ~isempty(session),
    try
        update_credentials(handles);
        status = 'connected';
    catch ME
        if isa(ME.ExceptionObject, 'Ice.CommunicatorDestroyedException')
            status = 'lost connection';
        else
            status = ME.message;
        end
    end
else
    status = 'not connected';
    
end
set(handles.text_status', 'String', sprintf('Status: %s', status));
set(handles.text_copyright, 'String', getLCCBCopyright())

% Choose default command line output for omeroLoginGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes omeroLoginGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = omeroLoginGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on key press with focus on edit_password and none of its controls.
function edit_password_KeyPressFcn(hObject, eventdata, handles)

password = get(hObject, 'UserData');
key = eventdata.Key;
switch key
    case 'backspace'
        password = password(1:end-1); % Delete the last character in the password
    case 'return'  % This cannot be done through callback without making tab to the same thing
        % do nothing
    case 'shift'
        % do nothing
    otherwise
        password = [password eventdata.Character];     
end

asterisk = password;
asterisk(1:end) = '*'; % Create a string of asterisks the same size as the password
set(hObject, 'String', asterisk,'UserData', password);
guidata(hObject,handles);


% --- Executes on button press in pushbutton_login.
function pushbutton_login_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_login (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

server  = get(handles.edit_server, 'String');
if isempty(server),
    errordlg('Please enter a valid server name', 'Server error', 'modal');
    return;
end

port  = get(handles.edit_port, 'String');
if isempty(port) || isnan(str2double(port)),
    errordlg('Please enter a valid port', 'Port error', 'modal');
    return;
end

username  = get(handles.edit_username, 'String');
if isempty(username),
    errordlg('Please enter a valid username', 'Username error', 'modal');
    return;
end

password  = get(handles.edit_password, 'UserData');
if isempty(password),
    errordlg('Please enter a valid password', 'Password error', 'modal');
    return;
end

% Create properties object to initialize the connection
properties = java.util.Properties();
properties.setProperty('omero.host', server);
properties.setProperty('omero.user', username);
properties.setProperty('omero.pass', password);
properties.setProperty('omero.port', port);
properties.setProperty('omero.keep_alive', '60');

connect(handles, properties);

% --- Executes on button press in pushbutton_connect_configuration_file.
function pushbutton_connect_configuration_file_Callback(hObject, eventdata, handles)


[file, path] = uigetfile('*.config', ['Select the configuration file to use to log in'...
    'to the OMERO server']);
if isequal(path, 0), return; end

connect(handles, [path file]);

function connect(handles, varargin)

global client
global session
try
    [client, session] = connectOmero(varargin{:});
    status = 'connected';
    update_credentials(handles)
catch ME
    if isa(ME.ExceptionObject, 'Ice.ConnectionRefusedException')
        status = 'connection refused';
    elseif isa(ME.ExceptionObject, 'Glacier2.PermissionDeniedException');
        status = 'password check failed';
    elseif isa(ME.ExceptionObject, 'Ice.DNSException');
        status = 'server name could not be resolved';
    else
        status = ME.message;
    end
end

set(handles.text_status, 'String', sprintf('Status: %s', status));


function update_credentials(handles)

global client
global session

adminService = session.getAdminService();
servername = char(client.getProperty('omero.host'));
set(handles.edit_server, 'String', servername)

userName = adminService.getEventContext().userName;
set(handles.edit_username, 'String', char(userName))