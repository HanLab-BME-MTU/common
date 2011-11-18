function varargout = correlationCalculationProcessGUI(varargin)
% correlationCalculationProcessGUI M-file for correlationCalculationProcessGUI.fig
%      correlationCalculationProcessGUI, by itself, creates a new correlationCalculationProcessGUI or raises the existing
%      singleton*.
%
%      H = correlationCalculationProcessGUI returns the handle to a new correlationCalculationProcessGUI or the handle to
%      the existing singleton*.
%
%      correlationCalculationProcessGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in correlationCalculationProcessGUI.M with the given input arguments.
%
%      correlationCalculationProcessGUI('Property','Value',...) creates a new correlationCalculationProcessGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before correlationCalculationProcessGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to correlationCalculationProcessGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help correlationCalculationProcessGUI

% Last Modified by GUIDE v2.5 18-Nov-2011 17:06:46

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @correlationCalculationProcessGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @correlationCalculationProcessGUI_OutputFcn, ...
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


% --- Executes just before correlationCalculationProcessGUI is made visible.
function correlationCalculationProcessGUI_OpeningFcn(hObject,eventdata,handles,varargin)

processGUI_OpeningFcn(hObject, eventdata, handles, varargin{:});

userData=get(handles.figure1,'UserData');
funParams = userData.crtProc.funParams_;

% Set up available input channels
set(handles.listbox_availableMovies,'String',userData.MD.movieDataFile_, ...
    'UserData',1:numel(userData.MD.movies_));

movieIndex = funParams.MovieIndex;

if ~isempty(movieIndex)
    movieString = userData.MD.movieDataFile_(movieIndex);
else
    movieString = {};
end

set(handles.listbox_selectedMovies,'String',movieString,...
    'UserData',movieIndex);

% Set up available input processes
allProcString = userData.crtProc.getCorrelationProcesses();
set(handles.listbox_availableProcesses,'String',allProcString);

% Set up selected input processes
procString = funParams.ProcessName;
set(handles.listbox_selectedProcesses,'String',procString);

set(handles.edit_BandMin,'String',funParams.BandMin);
set(handles.edit_BandMax,'String',funParams.BandMax);
userData.previewFig=-1;
% Choose default command line output for correlationCalculationProcessGUI
handles.output = hObject;

% Update user data and GUI data
set(hObject, 'UserData', userData);
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = correlationCalculationProcessGUI_OutputFcn(~, ~, handles) 
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

if isfield(userData, 'helpFig') && ishandle(userData.helpFig)
   delete(userData.helpFig) 
end

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

% --- Executes on button press in checkbox_allProcesses.
function checkbox_allProcesses_Callback(hObject, ~, handles)

availableProcesses = get(handles.listbox_availableProcesses, 'String');
if isempty(availableProcesses), return; end

if get(hObject,'Value')
    set(handles.listbox_selectedProcesses, 'String', availableProcesses);
else
    set(handles.listbox_selectedProcesses, 'String', {}, 'Value',1);
end

% --- Executes on button press in pushbutton_selectProcesses.
function pushbutton_selectProcesses_Callback(hObject, eventdata, handles)
% call back function of 'select' button

availableProcesses = get(handles.listbox_availableProcesses, 'String');
selectedProcesses = get(handles.listbox_selectedProcesses, 'String');
procID = get(handles.listbox_availableProcesses, 'Value');

newProcID = procID(~ismember(availableProcesses(procID),selectedProcesses));
selectedProcesses = vertcat(selectedProcesses,availableProcesses(newProcID));

set(handles.listbox_selectedProcesses, 'String', selectedProcesses);


% --- Executes on button press in pushbutton_deleteProcesses.
function pushbutton_deleteProcesses_Callback(hObject, eventdata, handles)
% Call back function of 'delete' button
selectedProcesses = get(handles.listbox_selectedProcesses,'String');
procID = get(handles.listbox_selectedProcesses,'Value');

% Return if list is empty
if isempty(selectedProcesses) || isempty(procID),return; end

% Refresh listbox
selectedProcesses(procID) = [ ];
set(handles.listbox_selectedProcesses,'String',selectedProcesses,...
    'Value',max(1,min(length(selectedProcesses),procID)));


% --- Executes on button press in checkbox_allMovies.
function checkbox_allMovies_Callback(hObject, eventdata, handles)


availableMovies = get(handles.listbox_availableMovies, 'String');
if isempty(availableMovies), return; end

availableMoviesIndex = get(handles.listbox_availableMovies, 'Userdata');

if get(hObject,'Value')
    set(handles.listbox_selectedMovies, 'String', availableMovies);
    selectedMoviesIndex = availableMoviesIndex;
else
    set(handles.listbox_selectedMovies, 'String', {}, 'Value',1);
    selectedMoviesIndex = [ ];
end
set(handles.listbox_selectedMovies, 'UserData', selectedMoviesIndex);


% --- Executes on button press in pushbutton_selectMovies.
function pushbutton_selectMovies_Callback(hObject, eventdata, handles)

% call back function of 'select' button
availableMovies = get(handles.listbox_availableMovies, 'String');
selectedMovies = get(handles.listbox_selectedMovies, 'String');
id = get(handles.listbox_availableMovies, 'Value');

% If channel has already been added, return;
availableMovieIndex = get(handles.listbox_availableMovies, 'Userdata');
selectedMovieIndex = get(handles.listbox_selectedMovies, 'Userdata');

for i = id
    if any(strcmp(availableMovies{i}, selectedMovies) )
        continue;
    else
        selectedMovies{end+1} = availableMovies{i};
        selectedMovieIndex = cat(2, selectedMovieIndex, availableMovieIndex(i));
    end
end

set(handles.listbox_selectedMovies, 'String', selectedMovies,...
    'Userdata', selectedMovieIndex);

% --- Executes on button press in pushbutton_deleteMovies.
function pushbutton_deleteMovies_Callback(hObject, eventdata, handles)
% Call back function of 'delete' button
selectedMovies = get(handles.listbox_selectedMovies,'String');
id = get(handles.listbox_selectedMovies,'Value');

% Return if list is empty
if isempty(selectedMovies) || isempty(id),return; end

% Delete selected item
selectedMovies(id) = [ ];

% Delete userdata
selectedMovieIndex = get(handles.listbox_selectedMovies, 'Userdata');
selectedMovieIndex(id) = [ ];
set(handles.listbox_selectedMovies, 'Userdata', selectedMovieIndex);

% Point 'Value' to the second last item in the list once the 
% last item has been deleted
if (id >length(selectedMovies) && id>1)
    set(handles.listbox_selectedMovies,'Value',length(selectedMovies));
end
% Refresh listbox
set(handles.listbox_selectedMovies,'String',selectedMovies);

% --- Executes on button press in checkbox_selectSlices.
function checkbox_selectSlices_Callback(hObject, eventdata, handles)

userData=get(handles.figure1,'UserData');
if ishandle(userData.previewFig), delete(userData.previewFig); end
if ~get(hObject,'Value'), return; end

movieProps = get(handles.listbox_selectedMovies,{'UserData','Value'});
movieID=movieProps{1}(movieProps{2});
p.ProcessName = get(handles.listbox_selectedProcesses,'String');
userData.previewFig=figure;

corrProc = CorrelationCalculationProcess(userData.MD.movies_{movieID},'');
parseProcessParams(corrProc,p);
input = corrProc.getInput;
nInput = numel(input);

alphamask = .2*ones(397,120);
alphamask(userData.crtProc.funParams_.SliceIndex{movieID},:)=1;
hAxes = -ones(nInput,1);
h = -ones(nInput,1);
for i=1:nInput;
    hAxes(i) = axes('Position',[(i-1)/nInput 0.1 1/nInput .8]);
    axesArgs={'hAxes',hAxes(i)};
    procID=input(i).processIndex;
    if ~isempty(input(i).channelIndex)
        chanArgs={input(i).channelIndex};
    else
        chanArgs={};
    end
    h(i) = userData.MD.movies_{movieID}.processes_{procID}.draw(chanArgs{:},axesArgs{:});
    hCbar = findobj(userData.previewFig,'Tag','Colorbar');
    delete(hCbar);    
    
    set(h(i),'AlphaData',alphamask,'AlphaDataMapping','none',...
        'ButtonDownFcn',@(h,event)editSlices(h,event,guidata(hObject)));
%     xLim = get(gca,'XLim');
%     r(i) = imrect(gca, [xLim(1) 10 xLim(2) 100]);
%     addNewPositionCallback(r(i),@(p) title(mat2str(p,3)));
%     fcn = makeConstrainToRectFcn('imrect',xLim,get(gca,'YLim'));
%     setPositionConstraintFcn(r(i),fcn);
%     
end
set(userData.previewFig,'DeleteFcn',@(h,event)closeGraphFigure(hObject));
set(handles.figure1, 'UserData', userData);

function editSlices(hObject)
set(hObject,'Value',0);   

 
function closeGraphFigure(hObject)
set(hObject,'Value',0);   


% --- Executes on button press in pushbutton_done.
function pushbutton_done_Callback(hObject, eventdata, handles)

% Check user input
if isempty(get(handles.listbox_selectedMovies, 'String'))
    errordlg('Please select at least one input process from ''Available Movies''.','Setting Error','modal')
    return;
end

if isempty(get(handles.listbox_selectedProcesses, 'String'))
    errordlg('Please select at least one input process from ''Available Processes''.','Setting Error','modal')
    return;
end

bandMin = str2double(get(handles.edit_BandMin,'String'));
if isnan(bandMin) || bandMin<1
    errordlg('Please enter a valid value for the layer of cells to be used',...
        'Setting error','modal');
    return;
end

bandMax = str2double(get(handles.edit_BandMax,'String'));
if isnan(bandMax) 
    errordlg('Please enter a valid value for the maximum layer of cells to be used',...
        'Setting error','modal');
    return;
end

funParams.MovieIndex = get(handles.listbox_selectedMovies, 'Userdata');
funParams.ProcessName = get(handles.listbox_selectedProcesses, 'String');
funParams.BandMin=bandMin;
funParams.BandMax=bandMax;

% Process Sanity check ( only check underlying data )
userData = get(handles.figure1, 'UserData');
try
    userData.crtProc.sanityCheck;
catch ME

    errordlg([ME.message 'Please double check your data.'],...
                'Setting Error','modal');
    return;
end

% Set parameters
processGUI_ApplyFcn(hObject, eventdata, handles,funParams);
