function varargout = movieDataGUI(varargin)
% MOVIEDATAGUI M-file for movieDataGUI.fig
%      MOVIEDATAGUI, by itself, creates a new MOVIEDATAGUI or raises the existing
%      singleton*.
%
%      H = MOVIEDATAGUI returns the handle to a new MOVIEDATAGUI or the handle to
%      the existing singleton*.
%
%      MOVIEDATAGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MOVIEDATAGUI.M with the given input arguments.
%
%      MOVIEDATAGUI('Property','Value',...) creates a new MOVIEDATAGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before movieDataGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to movieDataGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help movieDataGUI

% Last Modified by GUIDE v2.5 14-Nov-2011 13:42:46

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @movieDataGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @movieDataGUI_OutputFcn, ...
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


% --- Executes just before movieDataGUI is made visible.
function movieDataGUI_OpeningFcn(hObject, eventdata, handles, varargin)
%
% movieDataGUI('mainFig', handles.figure1) - call from movieSelector
% movieDataGUI(MD) - MovieData viewer
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
% NOTE: If movieDataGUI is under the "Overview" mode, additionally, 
% 
% userData.MD - the handle of selected MovieData object
%
%

% Input check
ip = inputParser;
ip.addRequired('hObject',@ishandle);
ip.addRequired('eventdata',@(x) isstruct(x) || isempty(x));
ip.addRequired('handles',@isstruct);
ip.addOptional('MD',[],@(x) isa(x,'MovieData'));
ip.addParamValue('mainFig',-1,@ishandle);
ip.parse(hObject,eventdata,handles,varargin{:})

% Store inpu
userData = get(handles.figure1, 'UserData');
userData.MD=ip.Results.MD;
userData.mainFig=ip.Results.mainFig;

[copyright openHelpFile] = userfcn_softwareConfig(handles);
set(handles.text_copyright, 'String', copyright)


% Set channel object array
userData.channels = [];

% Load help icon from dialogicons.mat
load lccbGuiIcons.mat
supermap(1,:) = get(hObject,'color');

userData.colormap = supermap;
userData.questIconData = questIconData;

set(handles.figure1,'CurrentAxes',handles.axes_help);
Img = image(questIconData);
set(hObject,'colormap',supermap);
set(gca, 'XLim',get(Img,'XData'),'YLim',get(Img,'YData'),...
    'visible','off');
set(Img,'ButtonDownFcn',@icon_ButtonDownFcn);

if openHelpFile
    set(Img, 'UserData', struct('class', mfilename))
end


if ~isempty(userData.MD),
    userData.channels = userData.MD.channels_;
        
    % Channel listbox
    cPath=arrayfun(@(x) x.channelPath_,userData.channels,'UniformOutput',false);
    set(handles.listbox_channel, 'String', cPath)
    
    % GUI setting
    set(handles.pushbutton_delete, 'Enable', 'off')
    set(handles.pushbutton_add, 'Enable', 'off')
    set(handles.pushbutton_output, 'Enable', 'off')
    
    set(hObject, 'Name', 'Movie Detail')
    set(handles.edit_path,'String', [userData.MD.movieDataPath_ filesep userData.MD.movieDataFileName_])
    set(handles.edit_output, 'String', userData.MD.outputDirectory_)
    set(handles.edit_notes, 'String', userData.MD.notes_)
    
    % GUI setting - parameters
    propNames={'pixelSize_','timeInterval_','numAperture_','camBitdepth_'};
    validProps = ~cellfun(@(x) isempty(userData.MD.(x)),propNames);
    
    propNames=propNames(validProps);
    cellfun(@(x) set(handles.(['edit_' x(1:end-1)]),'Enable','off',...
        'String',userData.MD.(x)),propNames)    
end

% Choose default command line output for movieDataGUI
handles.output = hObject;

% Update handles structure
set(handles.figure1,'UserData',userData)
guidata(hObject, handles);

% UIWAIT makes movieDataGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = movieDataGUI_OutputFcn(~, ~, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton_cancel.
function pushbutton_cancel_Callback(~, ~, handles)
% hObject    handle to pushbutton_cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
delete(handles.figure1);


% --- Executes on button press in pushbutton_done.
function pushbutton_done_Callback(~, ~, handles)

userData = get(handles.figure1,'UserData');

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
if isempty(outputDir) || ~exist(outputDir, 'dir')
    errordlg('Please provide a valid output path to save your results.', ...
               'Empty Output Path', 'modal');
    return;    
end

% Concatenate numerical parameters as movie options
propNames={'pixelSize_','timeInterval_','numAperture_','camBitdepth_'};
propHandles = cellfun(@(x) handles.(['edit_' x(1:end-1)]),propNames);
propStrings =get(propHandles,'String');
validProps = ~cellfun(@isempty,propStrings);
propNames=propNames(validProps);
propValues=num2cell(str2double(propStrings(validProps)))';
 
movieOptions = vertcat(propNames,propValues);
movieOptions = reshape(movieOptions,1,numel(propNames)*2);

% If movieDataGUI is under "Overview" mode
if ~isempty(get(handles.edit_notes, 'String'))
    movieOptions=horzcat(movieOptions,'notes_',get(handles.edit_notes, 'String'));
end

if ~isempty(userData.MD);
    % Overview mode - edit existing MovieDat
    try
        set(userData.MD,movieOptions{:});
    catch ME
        errormsg = sprintf([ME.message '.\n\Editing movie data failed.']);
        errordlg(errormsg, 'User Input Error','modal');
        return;
    end
    % Create a pointer to the MovieData object (to use the same
    % sanityCheck command later)
    MD=userData.MD; 
else
    % Create Movie Data
    try
        MD = MovieData(userData.channels, outputDir, movieOptions{:});
    catch ME
        errormsg = sprintf([ME.message '.\n\nCreating movie data failed.']);
        errordlg(errormsg, 'User Input Error','modal');
        return;
    end
end

try
    MD.sanityCheck
catch ME
    delete(MD);
    errormsg = sprintf('%s.\n\nPlease check your movie data. Movie data is not saved.',ME.message);
    errordlg(errormsg,'Channel Error','modal');
    return;
end

% Run the save method (should launch the dialog box asking for the object 
% path and filename)
MD.save(); 

% If new MovieData was created (from movieSelectorGUI)
if ishandle(userData.mainFig), 
    % Retrieve main window userData
    userData_main = get(userData.mainFig, 'UserData');
    
    % Check if files in movie list are saved in the same file
    handles_main = guidata(userData.mainFig);
    contentlist = get(handles_main.listbox_movie, 'String');
    movieDataFullPath = [MD.movieDataPath_ filesep MD.movieDataFileName_];
    if any(strcmp(movieDataFullPath, contentlist))
        errordlg('Cannot overwrite a movie data file which is already in the movie list. Please choose another file name or another path.','Error','modal');
        return
    end
    
    % Append  MovieData object to movie selector panel
    userData_main.MD = cat(2, userData_main.MD, MD);
    
    % Refresh movie list box in movie selector panel
    contentlist{end+1} = movieDataFullPath;
    nMovies = length(contentlist);
    set(handles_main.listbox_movie, 'String', contentlist, 'Value', nMovies)
    title = sprintf('Movie List: %s/%s movie(s)', num2str(nMovies), num2str(nMovies));
    set(handles_main.text_movie_1, 'String', title)
    
    % Save the main window data
    set(userData.mainFig, 'UserData', userData_main)
end
% Delete current window
delete(handles.figure1)


function edit_property_Callback(hObject, eventdata, handles)
% hObject    handle to edit_timeInterval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_timeInterval as text
%        str2double(get(hObject,'String')) returns contents of edit_timeInterval as a double

if ~isempty(get(hObject,'String'))
    propTag = get(hObject,'Tag');
    propName = [propTag(length('edit_')+1:end) '_'];
    propValue = str2double(get(hObject,'String'));
    if ~MovieData.checkValue(propName,propValue)
        warndlg('Invalid property value','Setting Error','modal');
        set(hObject,'BackgroundColor',[1 .8 .8]);
        return
    end
end
set(hObject,'BackgroundColor',[1 1 1])


% --- Executes on button press in pushbutton_delete.
function pushbutton_delete_Callback(hObject, eventdata, handles)

userData = get(handles.figure1, 'Userdata');

contents = get(handles.listbox_channel,'String');
% Return if list is empty
if isempty(contents), return; end
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
if ishandle(userData.mainFig), 
    handles_main = guidata(userData.mainFig);
    userData_main = get(handles_main.figure1, 'UserData');
    userDir =userData_main.userDir;
else
    userDir=pwd;
end

path = uigetdir(userDir, 'Add Channels ...');
if path == 0, return; end

% Get current list
contents = get(handles.listbox_channel,'String');
if any(strcmp(contents,path))
   warndlg('This directory has been selected! Please select a differenct directory.',...
       'Warning','modal');
   return; 
end

% Create path object and save it to userData
try
    newChannel= Channel(path);
    newChannel.sanityCheck();
catch ME
    errormsg = sprintf('%s.\n\nPlease check this is valid channel.',ME.message);
    errordlg(errormsg,'Channel Error','modal');
    return
end

userData.channels = cat(2, userData.channels, newChannel);
% Refresh listbox_channel
contents{end+1} = path;
set(handles.listbox_channel,'string',contents);

% Set user directory
sepDir = regexp(path, filesep, 'split');
dir = sepDir{1};
for i = 2: length(sepDir)-1
    dir = [dir filesep sepDir{i}];
end

if ishandle(userData.mainFig), 
    userData_main.userDir = dir;
    set(handles_main.figure1, 'UserData', userData_main)
end

set(handles.figure1, 'Userdata', userData)
guidata(hObject, handles);


% --- Executes during object deletion, before destroying properties.
function figure1_DeleteFcn(hObject, eventdata, handles)

userData = get(handles.figure1, 'UserData');

if isfield(userData, 'iconHelpFig') && ishandle(userData.iconHelpFig)
   delete(userData.iconHelpFig) 
end


% --- Executes on button press in pushbutton_output.
function pushbutton_output_Callback(hObject, eventdata, handles)

pathname = uigetdir(pwd,'Select a directory to store the processes output');
if isnumeric(pathname), return; end

set(handles.edit_output, 'String', pathname);


% --- Executes on button press in pushbutton_setting_chan.
function pushbutton_setting_chan_Callback(hObject, eventdata, handles)

userData = get(handles.figure1, 'UserData');

if isempty(userData.channels), return; end
assert( isa(userData.channels(1), 'Channel'), 'User-defined: Not a valid ''Channel'' object');

userData.setChannelFig = channelGUI('mainFig', handles.figure1, 'modal');

set(handles.figure1,'UserData',userData);
