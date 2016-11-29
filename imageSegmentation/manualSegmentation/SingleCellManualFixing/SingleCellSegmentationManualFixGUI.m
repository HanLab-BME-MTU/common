function varargout = SingleCellSegmentationManualFixGUI(varargin)
% SINGLECELLSEGMENTATIONMANUALFIXGUI MATLAB code for SingleCellSegmentationManualFixGUI.fig
%      SINGLECELLSEGMENTATIONMANUALFIXGUI, by itself, creates a new SINGLECELLSEGMENTATIONMANUALFIXGUI or raises the existing
%      singleton*.
%
%      H = SINGLECELLSEGMENTATIONMANUALFIXGUI returns the handle to a new SINGLECELLSEGMENTATIONMANUALFIXGUI or the handle to
%      the existing singleton*.
%
%      SINGLECELLSEGMENTATIONMANUALFIXGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SINGLECELLSEGMENTATIONMANUALFIXGUI.M with the given input arguments.
%
%      SINGLECELLSEGMENTATIONMANUALFIXGUI('Property','Value',...) creates a new SINGLECELLSEGMENTATIONMANUALFIXGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SingleCellSegmentationManualFixGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SingleCellSegmentationManualFixGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SingleCellSegmentationManualFixGUI

% Last Modified by GUIDE v2.5 14-Jan-2014 15:19:50

% this code was made by Hunter, adopted by Liya for single cell
% segmentation for sequences,2013

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SingleCellSegmentationManualFixGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @SingleCellSegmentationManualFixGUI_OutputFcn, ...
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


% --- Executes just before SingleCellSegmentationManualFixGUI is made visible.
function SingleCellSegmentationManualFixGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SingleCellSegmentationManualFixGUI (see VARARGIN)

% Choose default command line output for SingleCellSegmentationManualFixGUI
handles.output = hObject;

handles.segThisCell = [];
set(handles.chooseCell,'Enable','off')
set(handles.popupmenu_target_channel,'Enable','off')
set(handles.goForSegmentation,'Enable','off')
set(handles.CheckForSegmentation,'Enable','off')
set(handles.popupmenu_supplementary_channel,'Enable','off')
set(handles.popupmenu_single_cell_this_movie_ID,'Enable','off')


handles.target_channelIdx = 1;
handles.sup_channelIdx = 2;
handles.currCell   = 1;
% handles.currSingleCell   = 1;
handles.single_cell_ID = 1 ;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SingleCellSegmentationManualFixGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = SingleCellSegmentationManualFixGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in loadObject.
function loadObject_Callback(hObject, eventdata, handles)
% hObject    handle to loadObject (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[fileName,pathName] = uigetfile;
cellName{1} = 'Choose a movie';
set(handles.chooseCell,'String',cellName)
if ~isempty(fileName)
    mObj = load([pathName filesep fileName]);
    
    if isfield(mObj,'ML')
        
        ML = MovieList.load([pathName filesep fileName],0);
        
        nCell = numel(ML.movieDataFile_);
        for iCell = 1:nCell
            idx                      = max(regexp(ML.movies_{iCell}.movieDataPath_,filesep));
            cellName{iCell+1}         = ML.movies_{iCell}.movieDataPath_(idx+1:end);
        end
        handles.ML          = ML.movies_;
    elseif isfield(mObj,'MD')
        MD            = MovieData.load([pathName filesep fileName],0);
        nCell         = 1;
        idx           = max(regexp(MD.movieDataPath_,filesep));
        cellName{2}   = MD.movieDataPath_(idx+1:end);
        handles.ML{1} = MD;
    else
        error('This is not a movieObjet')
    end
    
    set(handles.chooseCell,'Enable','on')
    set(handles.goForSegmentation,'Enable','on')
    set(handles.CheckForSegmentation,'Enable','on')
   
    set(handles.chooseCell,'String',cellName)
    
    
    handles.cellName    = cellName;
    
    handles.loadedCells = false(nCell,1);
    handles.loadedMask  = cell(nCell,1);
    handles.loadedImage = cell(nCell,1);
    guidata(hObject, handles);
    
end

% --- Executes on selection change in chooseCell.
function chooseCell_Callback(hObject, eventdata, handles)
% hObject    handle to chooseCell (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns chooseCell contents as cell array
%        contents{get(hObject,'Value')} returns selected item from chooseCell
contents = cellstr(get(hObject,'String'));

handles.segThisCell = find(cell2mat(cellfun(@(x) strcmp(x,contents{get(hObject,'Value')}),handles.cellName,'Unif',0)));

if handles.segThisCell > 1
    set(handles.popupmenu_target_channel,'Enable','on');
    handles.currCell = handles.segThisCell - 1;
    
    currObj    = handles.ML{handles.currCell};
    nChannel   = numel(currObj.channels_);
    all_channelIdx = arrayfun(@num2str,1:nChannel,'Unif',0);
        
    set(handles.popupmenu_supplementary_channel,'Enable','on');
    set(handles.popupmenu_single_cell_this_movie_ID,'Enable','on');
    set(handles.CheckForSegmentation,'Enable','on');
    
    set(handles.popupmenu_target_channel,'String',all_channelIdx);
    set(handles.popupmenu_target_channel,'Value',1);
    set(handles.popupmenu_supplementary_channel,'String',all_channelIdx);
    
    handles.sup_channelIdx = min(2,nChannel);
    set(handles.popupmenu_supplementary_channel,'Value',handles.sup_channelIdx);
      
    guidata(hObject, handles);
end

% --- Executes during object creation, after setting all properties.
function chooseCell_CreateFcn(hObject, eventdata, handles)
% hObject    handle to chooseCell (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in goForSegmentation.
function goForSegmentation_Callback(hObject, eventdata, handles)
% hObject    handle to goForSegmentation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.segThisCell > 1
    currObj  = handles.ML{handles.currCell};
    currChan = handles.target_channelIdx;
    supChan = handles.sup_channelIdx;
    currSingleCell = handles.single_cell_ID;
    nFrame = currObj.nFrames_;
    
    if ~handles.loadedCells(handles.currCell)
        images = currObj.channels_(handles.target_channelIdx).loadImage(1:nFrame);
        maskIdx   = currObj.getProcessIndex('MaskRefinementProcess');
        cellMask  = arrayfun(@(x) currObj.processes_{maskIdx}.loadChannelOutput(currChan,x),1:nFrame,'Unif',0);
        sup_cellMask  = arrayfun(@(x) currObj.processes_{maskIdx}.loadChannelOutput(supChan,x),1:nFrame,'Unif',0);
        masks     = cat(3,cellMask{:}); 
        sup_masks     = cat(3,sup_cellMask{:}); 
        maskFileNames = currObj.processes_{maskIdx}.getOutMaskFileNames(currChan);
        orgPath = currObj.processes_{maskIdx}.outFilePaths_{currChan};
        handles.loadedCells(handles.currCell)  =  true;
        handles.loadedMask{handles.currCell}   = masks;
        handles.loadedSupMask{handles.currCell}   = sup_masks;
        handles.loadedImage{handles.currCell}  = images;
        handles.maskFileNames{handles.currCell}   = maskFileNames;
        handles.maskPath{handles.currCell}  = orgPath;
        guidata(hObject, handles);
    else
        masks  = handles.loadedMask{handles.currCell};
        images = handles.loadedImage{handles.currCell};
        sup_masks  = handles.loadedSupMask{handles.currCell};      
        try
            maskFileNames = handles.maskFileNames{handles.currCell};
            orgPath = handles.maskPath{handles.currCell};
        catch
            maskFileNames = [];
            orgPath = [];
        end
    end
  
    segIdx    = currObj.getPackageIndex('SegmentationPackage');
    
    if(isempty(segIdx))
         segIdx = currObj.getPackageIndex('FilamentAnalysisPackage');    
    end

    if(isempty(segIdx))
         segIdx = currObj.getPackageIndex('QFSMPackage');    
    end

    segPath   = currObj.packages_{segIdx}.outputDirectory_;
    
    
    truthPath = [segPath filesep 'FixedChannel' num2str(currChan) 'Cell' num2str(currSingleCell)];
    fixedPath = [segPath filesep 'OnlyFixedChannel' num2str(currChan) 'Cell' num2str(currSingleCell)];    
    compPath  = [segPath filesep 'completedFramesChannel' num2str(currChan) 'Cell' num2str(currSingleCell)];
    
    boxall = nan(size(masks,3),4);

    boxall(:,1)=1;
    boxall(:,3)=size(masks,1);
    boxall(:,2)=1;
    boxall(:,4)=size(masks,2);
    
    aux.isCompleted = [];
    
    if exist([compPath filesep 'completedFrames.mat'],'file')
        aux      = load([compPath filesep 'completedFrames.mat']);
        try
            boxall = aux.boxall;
        end
        
        goodMask = imDir([truthPath]);
        for iMask = find(aux.isCompleted)'
            try
                masks(:,:,iMask) = imread([truthPath filesep goodMask(iMask).name]);
            catch
                masks(:,:,iMask) = 0;
            end
        end
            
    end
    
    [outMasks,boxall,isCompleted] = manualSegmentationFixSingleCellTweakGUI(images,masks,sup_masks,[],aux.isCompleted,boxall,[],fixedPath,compPath,orgPath,maskFileNames{1});
    

%     iCell=currentSingleCellID;
%     truthPath = [segPath filesep 'FixedChannel' num2str(currChan) 'Cell' num2str(iCell)];
%     compPath  = [segPath filesep 'completedFramesChannel' num2str(currChan) 'Cell' num2str(iCell)];


    % do the complete checking first
    if ~exist(compPath,'dir')        
        try
            mkdir(compPath)        
        catch
            system(['mkdir -p ' compPath]);        
        end
    end
    
    save([compPath filesep 'completedFrames.mat'],'isCompleted','boxall');


    % then save the results frame by frame
    if ~exist(truthPath,'dir')        
        mkdir(truthPath)        
    end
    
    if ~exist(fixedPath,'dir')        
        mkdir(fixedPath)        
    end
    
    for iFrame = 1:nFrame
        imwrite(outMasks(:,:,iFrame),[truthPath filesep 'mask_' num2str(iFrame) '.tif'],'tif');
    end
    
    for iFrame = 1:nFrame
        if(isCompleted(iFrame)>0)
            imwrite(outMasks(:,:,iFrame),[fixedPath filesep 'mask_' num2str(iFrame) '.tif'],'tif');
        else
            % if not finished, delete the old marking
            if(exist([fixedPath filesep 'mask_' num2str(iFrame) '.tif'],'file'))
                delete([fixedPath filesep 'mask_' num2str(iFrame) '.tif']);                
            end
        end      
        
        % save all of them in the "truth" folder
        imwrite(outMasks(:,:,iFrame),[truthPath filesep 'mask_' num2str(iFrame) '.tif'],'tif');       
    end     
    
 end


% --- Executes on selection change in popupmenu_target_channel.
function popupmenu_target_channel_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_target_channel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_target_channel contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_target_channel
contents           = cellstr(get(hObject,'String'));
handles.target_channelIdx = str2num( contents{get(hObject,'Value')} );
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function popupmenu_target_channel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_target_channel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_single_cell_this_movie_ID.
function popupmenu_single_cell_this_movie_ID_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_single_cell_this_movie_ID (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_single_cell_this_movie_ID contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_single_cell_this_movie_ID

handles.single_cell_ID = get(hObject,'Value');
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function popupmenu_single_cell_this_movie_ID_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_single_cell_this_movie_ID (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
 set(hObject,'String',{'1','2','3','4','5'});
 set(hObject,'Value',1);
 


% --- Executes on selection change in popupmenu_supplementary_channel.
function popupmenu_supplementary_channel_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_supplementary_channel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_supplementary_channel contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_supplementary_channel
contents           = cellstr(get(hObject,'String'));
handles.sup_channelIdx = str2num( contents{get(hObject,'Value')} );
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function popupmenu_supplementary_channel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_supplementary_channel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in CheckForSegmentation.
function CheckForSegmentation_Callback(hObject, eventdata, handles)
% hObject    handle to CheckForSegmentation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.segThisCell > 1
    currObj  = handles.ML{handles.currCell};
    currChan = handles.target_channelIdx;
    supChan = handles.sup_channelIdx;
    currSingleCell = handles.single_cell_ID;
    nFrame = currObj.nFrames_;

    segIdx    = currObj.getPackageIndex('SegmentationPackage');
    if(isempty(segIdx))
        segIdx = currObj.getPackageIndex('FilamentAnalysisPackage');
    end
    segPath   = currObj.packages_{segIdx}.outputDirectory_;
    
%     
%     segPath(segPath=='\') = '/';
%     ind_L = segPath(segPath=='/');
    
    
    truthPath = [segPath filesep 'FixedChannel' num2str(currChan) 'Cell' num2str(currSingleCell)];
    fixedPath = [segPath filesep 'OnlyFixedChannel' num2str(currChan) 'Cell' num2str(currSingleCell)];    
    compPath  = [segPath filesep 'completedFramesChannel' num2str(currChan) 'Cell' num2str(currSingleCell)];
    
    
    if ~exist(fixedPath,'dir')
        if ~exist(truthPath,'dir')
            msgbox('There is no result for this movie this cell.')
        else
            
            if exist(compPath,'dir')
                load([compPath filesep 'completedFrames.mat'],'isCompleted','boxall');
                ind = find(isCompleted>0);
                
                mkdir(fixedPath);
                
%                 for iFrame = (ind(:))'                                   
%                     try
%                         copyfile([truthPath,filesep,'mask_',num2str(iFrame),'tif'],[fixedPath,filesep,'mask_',num2str(iFrame),'tif']);                    
%                     end
%                 end
%                 
%                 winopen(fixedPath); 
                winopen(truthPath);
            else                
                winopen(truthPath);
                msgbox('There is only old version of marking, which could include unfinished markings. Please save during marking to save all markings.')
            end
            
        end
    else
        if ispc
            winopen(fixedPath);
        else
            disp('Data are stored in: ');
            disp(fixedPath);
        end
    end    
    
    
    if exist(compPath,'dir')                        
        load([compPath filesep 'completedFrames.mat'],'isCompleted','boxall');
        ind = find(isCompleted>0);
        msgbox(['Completed frames for this movie this cell are: ', num2str(((ind(:))'))],'Previous Manual Fixing');
    end
    
 end
