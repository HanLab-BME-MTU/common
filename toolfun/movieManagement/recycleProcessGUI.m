function varargout = recycleProcessGUI(varargin)
% RECYCLEPROCESSGUI M-file for recycleProcessGUI.fig
%      RECYCLEPROCESSGUI, by itself, creates a new RECYCLEPROCESSGUI or raises the existing
%      singleton*.
%
%      H = RECYCLEPROCESSGUI returns the handle to a new RECYCLEPROCESSGUI or the handle to
%      the existing singleton*.
%
%      RECYCLEPROCESSGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in RECYCLEPROCESSGUI.M with the given input arguments.
%
%      RECYCLEPROCESSGUI('Property','Value',...) creates a new RECYCLEPROCESSGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before recycleProcessGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to recycleProcessGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help recycleProcessGUI

% Last Modified by GUIDE v2.5 24-Sep-2010 13:47:24

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @recycleProcessGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @recycleProcessGUI_OutputFcn, ...
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


% --- Executes just before recycleProcessGUI is made visible.
function recycleProcessGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% 
% recycleProcessGUI(process, package, 'mainFig', handles.figure1)
% 
% Input:
%
%   process - the array of processes for recycle
%   package - the package where the processes would be attached to
%
% User Data:
% 
% userData.process - the array of processes for recycle
% userData.pacakge - the package where the processes would be attached to
% 
% userData.mainFig - handle of movie selector GUI
% 
% 

set(handles.text_copyright, 'String', userfcn_copyright)

userData = get(handles.figure1, 'UserData');
% Choose default command line output for recycleProcessGUI
handles.output = hObject;


if nargin > 3
    
    assert( all( cellfun(@(x)isa(x, 'Process'), varargin{1}) ), 'User-defined: The first input must be a cell array containing Process objects.')
    userData.process = varargin{1};
    
    assert(isa(varargin{2}, 'Package'), 'User-defined: The second input must be a Package object where the processes would be attached to.')
    userData.package = varargin{2};
    
    t = find(strcmp(varargin, 'mainFig'));
    assert( ~isempty(t), 'User-defined: Need to pass the handle of main figure as input.')
    userData.mainFig = varargin{t+1};

    
    % GUI set-up
    set(handles.text_package, 'String', userData.package.name_)
    set(handles.text_movie, 'String', [userData.package.owner_.movieDataPath_ userData.package.owner_.movieDataFileName_])
    
    string = cell(1, length(userData.process));
    for i = 1: length(userData.process)
        string{i} = [96+i '. ' userData.process{i}.name_ ' Step'];
    end
    set(handles.listbox_1, 'String', string, 'UserData', userData.process)
    
else
    error('User-defined: Not enough input arguments.')
end


set(handles.figure1,'UserData',userData)

uiwait(handles.figure1)
guidata(hObject, handles);

% UIWAIT makes recycleProcessGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = recycleProcessGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

delete(handles.figure1)


% --- Executes on button press in pushbutton_cancel.
function pushbutton_cancel_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
userData = get(handles.figure1, 'UserData');

uiresume(handles.figure1)


% --- Executes on button press in pushbutton_done.
function pushbutton_done_Callback(hObject, eventdata, handles)

% Add the processes to the new package
userData = get(handles.figure1, 'UserData');

userData_list2 = get(handles.listbox_2, 'UserData');

if isempty(userData_list2)
    errordlg('No step is selected.', 'modal')
    return 
end

processClassNames = userData.package.processClassNames_;

check = [];

for i = 1: length(userData_list2)
       
    id = find(cellfun(@(x)isa(userData_list2{i}, x),processClassNames ));
    
    if length(id) > 1
        % Special cases: Mask Refinement Process
        if isa(userData_list2{i}, 'MaskRefinementProcess')
            
            id = find(cellfun(@(x)strcmp(x, 'MaskRefinementProcess'),processClassNames ));
            userData.package.setProcess(id, userData_list2{i})
        else
            error('User-defined: process belongs to more than 1 classes in package''s process list.')
        end
        
    elseif length(id) ==1
        % General cases: child class, specific class
        userData.package.setProcess(id, userData_list2{i})
    else
        error('User-defined: Inapproprate process class in listbox2.')
    end
    
    % At last, double check to make sure package.processes_{id} is set no
    % more than once
    check = cat(2, check, id);
    if length(check) ~= length(unique(check))
        error('User-defined: one of the process list value in current package is changed more than once in this session.')
    end
end


uiresume(handles.figure1)

% --- Executes on selection change in listbox_1.
function listbox_1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_1


% --- Executes during object creation, after setting all properties.
function listbox_1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_add.
function pushbutton_add_Callback(hObject, eventdata, handles)


contentlist1 = get(handles.listbox_1, 'String');
contentlist2 = get(handles.listbox_2, 'String');

if isempty(contentlist1)
    return
end

userData_list1 = get(handles.listbox_1, 'UserData');
userData_list2 = get(handles.listbox_2, 'UserData');

id = get(handles.listbox_1, 'value');

if ~isempty(userData_list2)

% Avoid adding the same class of objects
if isa(userData_list1{id}, 'SegmentationProcess') && ~isa(userData_list1{1}, 'MaskRefinementProcess')

    if any(xor (cellfun(@(x)isa(x, 'SegmentationProcess'), userData_list2), ...
            cellfun(@(x)isa(x, 'MaskRefinementProcess'), userData_list2)))
       
        errordlg('Segmentation step has already been selected.', 'modal')
        return
    end    
else
    
    temp = class(userData_list1{id});
    if any(cellfun(@(x)isa(x, temp), userData_list2) )
       
        errordlg('This step has already been selected.', 'modal')
        return
    end
    
end

end

contentlist2{end+1} = contentlist1{id};
userData_list2{end +1} = userData_list1{id};

contentlist1(id) = [];
userData_list1(id) = [];

if (id > length(contentlist1) && id > 1)
    set(handles.listbox_1, 'Value', length(contentlist1));
end

set(handles.listbox_1, 'String', contentlist1, 'UserData', userData_list1)
set(handles.listbox_2, 'String', contentlist2, 'UserData', userData_list2)



% --- Executes on selection change in listbox_2.
function listbox_2_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_2 contents as cell array
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


% --- Executes on button press in pushbutton_remove.
function pushbutton_remove_Callback(hObject, eventdata, handles)

contentlist1 = get(handles.listbox_1, 'String');
contentlist2 = get(handles.listbox_2, 'String');

if isempty(contentlist2)
    return
end

userData_list1 = get(handles.listbox_1, 'UserData');
userData_list2 = get(handles.listbox_2, 'UserData');

id = get(handles.listbox_2, 'value');

contentlist1{end + 1} = contentlist2{id};
userData_list1{end+1} = userData_list2{id};

contentlist2(id) = [];
userData_list2(id) = [];

if (id > length(contentlist2) && id > 1)
    set(handles.listbox_2, 'Value', length(contentlist2));
end

set(handles.listbox_1, 'String', contentlist1, 'UserData', userData_list1)
set(handles.listbox_2, 'String', contentlist2, 'UserData', userData_list2)


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isequal(get(handles.figure1, 'waitstatus'), 'waiting')
    % The GUI is still in UIWAIT, us UIRESUME
    uiresume(handles.figure1);
else
    % The GUI is no longer waiting, just close it
    delete(handles.figure1);
end


% --- Executes on button press in pushbutton_detail.
function pushbutton_detail_Callback(hObject, eventdata, handles)

if isempty(get(handles.listbox_1, 'String'))
    return
end
userData = get(handles.figure1, 'UserData');

id = get(handles.listbox_1, 'Value');
userData_list1 = get(handles.listbox_1, 'UserData');
process = userData_list1{id};

% if newMovieDataGUI exist
if isfield(userData, 'detailFig') && ishandle(userData.detailFig)
    delete(userData.detailFig)
end

switch class(process)
    
    case 'ThresholdProcess'
        
        userData.detailFig = viewThresholdProcessGUI(process, 'mainFig', handles.figure1);
        
    case 'MaskRefinementProcess'
        
        userData.detailFig = viewMaskRefinementProcessGUI(process, 'mainFig', handles.figure1);
        
    otherwise
        
        msg = sprintf('The selected % step cannot be previewed at this time.', class(process.name_));
        errordlg(msg, 'modal');
        return
end

set(handles.figure1,'UserData',userData);


% --- Executes during object deletion, before destroying properties.
function figure1_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
userData = get(handles.figure1, 'UserData');

if isfield(userData, 'detailFig') && ishandle(userData.detailFig)
   delete(userData.detailFig) 
end
