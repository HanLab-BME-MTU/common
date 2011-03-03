function varargout = movieSelectorGUI(varargin)
% MOVIESELECTORGUI M-file for movieSelectorGUI.fig
%      MOVIESELECTORGUI, by itself, creates a new MOVIESELECTORGUI or raises the existing
%      singleton*.
%
%      H = MOVIESELECTORGUI returns the handle to a new MOVIESELECTORGUI or the handle to
%      the existing singleton*.
%
%      MOVIESELECTORGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MOVIESELECTORGUI.M with the given input arguments.
%
%      MOVIESELECTORGUI('Property','Value',...) creates a new MOVIESELECTORGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before movieSelectorGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to movieSelectorGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help movieSelectorGUI

% Last Modified by GUIDE v2.5 15-Aug-2010 13:34:44

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @movieSelectorGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @movieSelectorGUI_OutputFcn, ...
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


% --- Executes just before movieSelectorGUI is made visible.
function movieSelectorGUI_OpeningFcn(hObject, eventdata, handles, varargin)
%
% Useful tools:
% 
% User Data:
%
%   userData.MD - new or loaded MovieData object
%   userData.ML - newly saved or loaded MovieList object
%
%   userData.userDir - default open directory
%   userData.colormap - color map (used for )
%   userData.questIconData - image data of question icon
%
%   userData.firstPackageGUI - the name of package GUI
%
%   userData.newFig - handle of new movie set-up GUI
%   userData.iconHelpFig - handle of help dialog
%   userData.msgboxGUI - handle of message box GUI
%   userData.relocateFig - handle of relocateMovieDataGUI
%

[copyright openHelpFile] = userfcn_softwareConfig(handles);
set(handles.text_copyright, 'String', copyright)

userData = get(handles.figure1, 'UserData');

% Choose default command line output for setupMovieDataGUI
handles.output = hObject;


% set(handles.radiobutton_package_3, 'Visible', 'off')  %%%%%%%%%%%%%%%%%%


% Set callback function of radio button group uipanel_package
set(handles.uipanel_package, 'SelectionChangeFcn', @uipanel_package_SelectionChangeFcn);
 
% other user data set-up
userData.MD = [ ];
userData.ML = [ ];
userData.userDir = pwd;

% Indicate the first package control panel the GUI will go to after 
% MovieData is sucessfully created and MovieDate setup panel is destroied.
userData.firstPackageGUI = { @segmentationPackageGUI, ...
                             @biosensorsPackageGUI, ...
                             @uTrackPackageGUI};
                                
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
    set(Img, 'UserData', struct('class', 'MovieData'))
else
    set(Img, 'UserData', MovieData.getHelp)
end

% Save userdata
set(handles.figure1,'UserData',userData);
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = movieSelectorGUI_OutputFcn(hObject, eventdata, handles) 
% %varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton_cancel.
function pushbutton_cancel_Callback(hObject, eventdata, handles)
delete(handles.figure1)


% --- Executes on button press in pushbutton_done.
function pushbutton_done_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_done (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isempty( get(handles.listbox_movie, 'String') )
   warndlg('Please select at least one movie to continue.', 'Movie Selector', 'modal')
   return
end

if isempty(get(handles.uipanel_package, 'SelectedObject'))
   warndlg('Please select a package to continue.', 'Movie Selector', 'modal')
   return
end

userData = get(handles.figure1, 'userdata');
switch get(get(handles.uipanel_package, 'SelectedObject'), 'tag')

    
    case 'radiobutton_package_1'
        userData.firstPackageGUI{1}(userData.MD)
        delete(handles.figure1)
        
    case 'radiobutton_package_2'
        userData.firstPackageGUI{2}(userData.MD)
        delete(handles.figure1)
        
    case 'radiobutton_package_3'
        userData.firstPackageGUI{3}(userData.MD)
        delete(handles.figure1)        
        
end

% --- Executes on selection change in listbox_movie.
function listbox_movie_Callback(hObject, eventdata, handles)

contentlist = get(handles.listbox_movie, 'String');

if isempty(contentlist)
    title = sprintf('Movie List: 0/0 movie(s)');
    set(handles.text_movie_1, 'String', title)
else
    title = sprintf('Movie List: %s/%s movie(s)', num2str(get(handles.listbox_movie, 'Value')), num2str(length(contentlist)));
    set(handles.text_movie_1, 'String', title)    
end

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_movie contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_movie


% --- Executes during object creation, after setting all properties.
function listbox_movie_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_movie (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in pushbutton_new.
function pushbutton_new_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_new (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
userData = get(handles.figure1, 'UserData');

% if newMovieDataGUI exist
if isfield(userData, 'newFig') && ishandle(userData.newFig)
    delete(userData.newFig)
end
userData.newFig = newMovieDataGUI('mainFig',handles.figure1);
set(handles.figure1,'UserData',userData);


% --- Executes on button press in pushbutton_delete.
function pushbutton_delete_Callback(hObject, eventdata, handles)

userData = get(handles.figure1, 'Userdata');

contentlist = get(handles.listbox_movie,'String');
% Return if list is empty
if isempty(contentlist)
    return;
end
num = get(handles.listbox_movie,'Value');

% Delete channel object
delete(userData.MD(num))
userData.MD(num) = [];

% Refresh listbox_channel
contentlist(num) = [ ];
set(handles.listbox_movie,'String',contentlist);

% Point 'Value' to the second last item in the list once the 
% last item has been deleted
if num>length(contentlist) && num>1
    set(handles.listbox_movie,'Value',length(contentlist));
end
if isempty(contentlist)
    title = sprintf('Movie List: 0/0 movie(s)');
    set(handles.text_movie_1, 'String', title)
else
    title = sprintf('Movie List: %s/%s movie(s)', num2str(get(handles.listbox_movie, 'Value')), num2str(length(contentlist)));
    set(handles.text_movie_1, 'String', title)    
end

set(handles.figure1, 'Userdata', userData)
guidata(hObject, handles);



% --- Executes on button press in pushbutton_detail.
function pushbutton_detail_Callback(hObject, eventdata, handles)
if isempty(get(handles.listbox_movie, 'String'))
    return
end
userData = get(handles.figure1, 'UserData');

% if newMovieDataGUI exist
if isfield(userData, 'newFig') && ishandle(userData.newFig)
    delete(userData.newFig)
end

userData.newFig = newMovieDataGUI('mainFig',handles.figure1, 'overview', userData.MD(get(handles.listbox_movie, 'value')));
set(handles.figure1,'UserData',userData);


function uipanel_package_SelectionChangeFcn(hObject, eventdata)

% function userfcn_display(handles, method)
% Function to enable/disable uicontrols

% userData = get(handles.figure1,'UserData');

% switch method
%     case 'method1'
        
        % Enable 'New Movie Data'
%         set(handles.pushbutton_new, 'Enable', 'on')
%          set(handles.pushbutton_detail, 'Enable', 'on')
%         set(handles.pushbutton_delete, 'Enable', 'on')
%         set(handles.listbox_movie, 'Enable', 'on', 'BackgroundColor','white');

        
        % Disable 'Existing Movie Data'
%         set(handles.pushbutton_open, 'Enable', 'off')
%         set(handles.pushbutton_detail_load, 'Enable', 'off')
%         set(handles.listbox_load, 'Enable', 'off', 'BackgroundColor',[.95 .95 .95]);  
        
%         set(handles.pushbutton_done,'string','Save & Apply');
       
        
%     case 'method2'
        
        % Enable 'Existing Movie Data'
%         set(handles.pushbutton_open, 'Enable', 'on')
%         set(handles.pushbutton_detail_load, 'Enable', 'on')
%         set(handles.listbox_load, 'Enable', 'on', 'BackgroundColor','white');
        
        % Disable 'New Movie Data'
%         set(handles.pushbutton_new, 'Enable', 'off')
%         set(handles.pushbutton_detail, 'Enable', 'off')
%         set(handles.pushbutton_delete, 'Enable', 'off')
%         set(handles.listbox_movie, 'Enable', 'off', 'BackgroundColor',[.95 .95 .95])
        
%         set(handles.pushbutton_done,'string','Apply');
       
        
%     otherwise
%         error('User-defined: Error using function userfcn_display')
% end
            


% --- Executes during object deletion, before destroying properties.
function figure1_DeleteFcn(hObject, eventdata, handles)
% Delete function

userData = get(handles.figure1, 'UserData');

% Delete new window
if isfield(userData, 'newFig')
        if userData.newFig~=0 && ishandle(userData.newFig)
            delete(userData.newFig)
        end
end

if isfield(userData, 'iconHelpFig') && ishandle(userData.iconHelpFig)
   delete(userData.iconHelpFig) 
end
if isfield(userData, 'msgboxGUI') && ishandle(userData.msgboxGUI)
   delete(userData.msgboxGUI)
end
if isfield(userData, 'relocateFig') && ishandle(userData.relocateFig)
   delete(userData.relocateFig)
end



% --- Executes on button press in pushbutton_open.
function pushbutton_open_Callback(hObject, eventdata, handles)

userData = get(handles.figure1, 'UserData');

[filename, pathname] = uigetfile('*.mat','Select Movie Data MAT file', userData.userDir);
if ~any([filename pathname])
    return;
else
    userData.userDir = pathname;
end

% Check if reselect the movie data that is already in the listbox
contentlist = get(handles.listbox_movie, 'String');

if any(strcmp([pathname filename], contentlist))
    errordlg('This movie data has already been selected.','Error','modal');
    return
end

try
    pre = whos('-file', [pathname filename]);  % - Exception: fail to access .mat file
catch ME
    errordlg(ME.message,'Fail to open MAT file.','modal');
    return;
end

% Find MovieData object in .mat file before loading it
structM1 = pre( logical(strcmp({pre(:).class},'MovieData')) );
structM2 = pre( logical(strcmp({pre(:).class},'MovieList')) );

if ~isempty(structM1)
    structM = structM1;
    type = 'MovieData';
    
elseif ~isempty(structM2)
    structM = structM2;
    type = 'MovieList';
    
else
    errordlg('No movie data or movie list is found in selected MAT-file.',...
            'MAT File Error','modal');
    return
end

if length(structM)>1
    
    % Exception - multiple movies saved in one MAT file
    errordlg('The selected MAT-file contains more than one movie data or movie list.',...
            'MAT File Error','modal');
    return
else
    % Load MovieData or MovieList object in .mat file
    load([pathname filename],'-mat',structM.name);
    % Name/Rename MovieData obj as 'MD'
    eval(['M =' structM.name ';']);
end

switch type
    case 'MovieData'
        
        try 
            M.sanityCheck(pathname, filename)
        catch ME
            
            if isfield(userData, 'newFig') && ishandle(userData.newFig)
                delete(userData.newFig)
            end
            userData.newFig = newMovieDataGUI('mainFig',handles.figure1, 'overview', M);
            msg = sprintf('Movie Data: %s\n\nError: %s\n\nMovie data is not successfully loaded. Please refer to movie detail and adjust your data.', [pathname filename],ME.message);
            errordlg(msg, 'Movie Data Error','modal'); 
            return
        end
        userData.MD = cat(2, userData.MD, M);
        contentlist{end+1} = [pathname filename];
        
    case 'MovieList'
        % Find duplicate movie data in list box
        movieDataFile = M.movieDataFile_;
        index = 1: length(movieDataFile);
        index = index( ~cellfun(@(z)any(z), cellfun(@(x)strcmp(x, contentlist), movieDataFile, 'UniformOutput', false), 'UniformOutput', true) );
        
        if isempty(index)
            msg = sprintf('All movie data in movie list file %s has already been added to the movie list box.', M.movieListFileName_);
            warndlg(msg,'Warning','modal');
            return
        end
        
        reloadME = [];
        errorME = [];
        healthMD = [];
        
        [movieException, MDList] = M.sanityCheck(index, pathname, filename );
        
        % Explore cell array 'movieException'
        for i = 1: length(movieException)
            
            if isempty(movieException{i}) 
                healthMD = cat(2, healthMD, i);
            elseif strcmp(movieException{i}.identifier, 'lccb:ml:nofile' );
                reloadME = cat(2, reloadME, i);
            else
                errorME = cat(2, errorME, i);
            end
        end
        
        % Error movie index
        if ~isempty(errorME)
            filemsg = '';
            for i = errorME
                filemsg = cat(2, filemsg, sprintf('Movie %d:  %s\nError:  %s\n\n', index(i), movieDataFile{index(i)}, movieException{i}.message));
            end
            msg = sprintf('The following movie(s) cannot be sucessfully loaded:\n\n%s', filemsg);
            titlemsg = sprintf('Movie List: %s', [pathname filename]);
            userData.msgboxGUI = msgboxGUI('title',titlemsg,'text', msg);
        end
        
        % Healthy Movie Data
        if ~isempty(healthMD)
            userData.ML = cat(2, userData.ML, M);
            for i = healthMD
                userData.MD = cat(2, userData.MD, MDList{i});
            end
            contentlist = cat(2, contentlist', movieDataFile(index(healthMD)) );
        end
        
        % Reload movie index
        if ~isempty(reloadME)
            filemsg = '';
            for i = reloadME
                filemsg = cat(2, filemsg, sprintf('\n%s', movieDataFile{index(i)} ) );
            end
            
            % Ask user if relocate the movie data files
            msg = sprintf('Cannot find file(s):\n%s\n\nDo you want to relocate the above files?', filemsg);
            user_response = questdlg(msg, 'Movie Listbox', 'Yes','No','Yes');
            
            switch lower(user_response)
                case 'yes'
                    userData.relocateFig = relocateMovieDataGUI(M, index(reloadME), 'mainFig', handles.figure1);
                case 'no'
                    msg = sprintf('Do you want to keep or remove the unfound movie(s) in movielist %s?\n%s', filename, filemsg);
                    user_response2 = questdlg(msg, 'Movie Listbox', 'Keep Movie(s)', 'Remove Movie(s)', 'Keep Movie(s)');
                    
                    if strcmpi( user_response2, 'remove movie(s)')

                            M.removeMovieDataFile(index(reloadME));
                            M.saveMovieList
                            msg = sprintf('The following movie(s) have been removed from movielist %s.\n%s', filename, filemsg);
                            msgbox(msg, 'Help', 'help', 'modal');
                    end
            end
            
        end        
        
    otherwise
        error('User-defined: varable ''type'' does not have an approprate value.')
end

% Refresh movie list box in movie selector panel

set(handles.listbox_movie, 'String', contentlist, 'Value', length(contentlist))
title = sprintf('Movie List: %s/%s movie(s)', num2str(get(handles.listbox_movie, 'Value')), num2str(length(contentlist)));
set(handles.text_movie_1, 'String', title)


set(handles.figure1, 'UserData', userData);


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
function menu_file_new_Callback(hObject, eventdata, handles)

pushbutton_new_Callback(handles.pushbutton_new, [], handles)


% --------------------------------------------------------------------
function menu_file_open_Callback(hObject, eventdata, handles)

pushbutton_open_Callback(handles.pushbutton_open, [], handles)


% --------------------------------------------------------------------
function menu_file_quit_Callback(hObject, eventdata, handles)
% hObject    handle to menu_file_quit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
delete(handles.figure1)


% --- Executes on button press in pushbutton_deleteall.
function pushbutton_deleteall_Callback(hObject, eventdata, handles)
userData = get(handles.figure1, 'Userdata');

contentlist = get(handles.listbox_movie,'String');
% Return if list is empty
if isempty(contentlist)
    return;
end
 
user_response = questdlg(['Are you sure to delete all the ' num2str(length(contentlist)) ' movie(s) in the listbox?'], ...
    'Movie Listbox', 'Yes','No','Yes');

if strcmpi('no', user_response)
    return
end

% Delete channel object
userData.MD = [];

% Refresh listbox_channel
contentlist = {};
set(handles.listbox_movie,'String',contentlist, 'Value',1);

    title = sprintf('Movie List: 0/0 movie(s)');
    set(handles.text_movie_1, 'String', title)

set(handles.figure1, 'Userdata', userData)
guidata(hObject, handles);


% --- Executes on button press in pushbutton_save.
function pushbutton_save_Callback(hObject, eventdata, handles)

userData = get(handles.figure1, 'UserData');

contentList = get(handles.listbox_movie, 'String');
if isempty(contentList)
    warndlg('No movie selected. Please create new movie data or open existing movie data or movie list.', 'No Movie Selected', 'modal')
    return
end

if isempty(userData.ML)
    movieListPath = [userData.userDir filesep];
    movieListFileName = 'movieDataList.mat';
else
    movieListPath = userData.ML(end).movieListPath_;
    movieListFileName = userData.ML(end).movieListFileName_;
end

% Ask user where to save the movie data file
[file,path] = uiputfile('*.mat','Find a place to save your movie data',...
             [movieListPath filesep movieListFileName]);
         
if ~any([file,path])
    return
end

try
    ML = MovieList(contentList, path, file);
catch ME
    msg = sprintf('%s\n\nMovie list is not saved.', ME.message);
    errordlg(msg, 'Movie List Error', 'modal')
    return
end

ML.saveMovieList
