function mainFig = movieDataViewer(varargin)

ip = inputParser;
ip.addRequired('MD',@(x) isa(x,'MovieData'));
ip.addOptional('procId',0,@isnumeric);
ip.parse(varargin{:});

% Chek
h=findobj(0,'Name','Movie Viewer');
if ~isempty(h), delete(h); end
mainFig=figure('Name','Movie Viewer','Position',[0 0 200 200],...
    'NumberTitle','off','Tag','figure1','Toolbar','none','MenuBar','none',...
    'Color',get(0,'defaultUicontrolBackgroundColor'),'Resize','off',...
    'DeleteFcn', @(h,event) deleteViewer(guidata(h)));
userData=get(mainFig,'UserData');
userData.MD=ip.Results.MD;

% Classify processes by type (image or overlay)
validProcId= find(cellfun(@(x) ismember('getDrawableOutput',methods(x)) &...
    x.success_,userData.MD.processes_));
validProc=userData.MD.processes_(validProcId);
isImageProc =cellfun(@(x) any(strcmp({x.getDrawableOutput.type},'image')),validProc);
imageProc=validProc(isImageProc);
imageProcId = validProcId(isImageProc);
isOverlayProc =cellfun(@(x) any(strcmp({x.getDrawableOutput.type},'overlay')),validProc);
overlayProc=validProc(isOverlayProc);
overlayProcId = validProcId(isOverlayProc);

% Create series of anonymous function to generate process controls
createProcText= @(panel,i,j,pos,name) uicontrol(panel,'Style','text',...
    'Position',[10 pos 200 20],'Tag',['text_process' num2str(i)],...
    'String',name,'HorizontalAlignment','left','FontWeight','bold');
createOutputText= @(panel,i,j,pos,text) uicontrol(panel,'Style','text',...
    'Position',[40 pos 200 20],'Tag',['text_process' num2str(i) '_output'...
    num2str(j)],'String',text,'HorizontalAlignment','left');
createChannelButton= @(panel,i,j,k,pos) uicontrol(panel,'Style','radio',...
    'Position',[200+30*k pos 20 20],'Tag',['radiobutton_process' num2str(i) '_output'...
    num2str(j) '_channel' num2str(k)]);
createChannelBox= @(panel,i,j,k,pos) uicontrol(panel,'Style','checkbox',...
    'Position',[200+30*k pos 20 20],'Tag',['checkbox_process' num2str(i) '_output'...
    num2str(j) '_channel' num2str(k)],...
    'Callback',@(h,event) redrawOverlay(h,guidata(h)));

%% Create image panel
imagePanel = uibuttongroup(mainFig,'Position',[0 0 1/2 1],...
    'Title','Image','BackgroundColor',get(0,'defaultUicontrolBackgroundColor'),...
    'Units','pixels','Tag','uipanel_image');

% Create image option controls
hPosition1=10;
if isempty(userData.MD.timeInterval_),
    timeStampStatus = 'off';
else
    timeStampStatus = 'on';
end
uicontrol(imagePanel,'Style','checkbox',...
    'Position',[10 hPosition1 200 20],'Tag','checkbox_timeStamp',...
    'String',' Time stamp','HorizontalAlignment','left','FontWeight','bold',...
    'Enable',timeStampStatus,...
    'Callback',@(h,event) setTimeStamp(guidata(h)));
uicontrol(imagePanel,'Style','popupmenu','Position',[150 hPosition1 120 20],...
    'String',{'NorthEast', 'SouthEast', 'SouthWest', 'NorthWest'},'Value',4,...
    'Tag','popupmenu_timeStampLocation',...
    'Enable',timeStampStatus,...
    'Callback',@(h,event) setTimeStamp(guidata(h)));

hPosition1=hPosition1+30;
if isempty(userData.MD.pixelSize_),
    scaleBarStatus = 'off';
else
    scaleBarStatus = 'on';
end
uicontrol(imagePanel,'Style','edit','Position',[40 hPosition1 50 20],...
    'String','1','BackgroundColor','white','Tag','edit_imageScaleBar',...
    'Enable',scaleBarStatus,...
    'Callback',@(h,event) setScaleBar(guidata(h),'imageScaleBar'));
uicontrol(imagePanel,'Style','text','Position',[100 hPosition1-2 70 20],...
    'String','microns','HorizontalAlignment','left');

uicontrol(imagePanel,'Style','checkbox',...
    'Position',[180 hPosition1 100 20],'Tag','checkbox_imageScaleBarLabel',...
    'String',' Show label','HorizontalAlignment','left',...
    'Enable',scaleBarStatus,...
    'Callback',@(h,event) setScaleBar(guidata(h),'imageScaleBar'));

hPosition1=hPosition1+30;
uicontrol(imagePanel,'Style','checkbox',...
    'Position',[10 hPosition1 200 20],'Tag','checkbox_imageScaleBar',...
    'String',' Scalebar','HorizontalAlignment','left','FontWeight','bold',...
    'Enable',scaleBarStatus,...
    'Callback',@(h,event) setScaleBar(guidata(h),'imageScaleBar'));
uicontrol(imagePanel,'Style','popupmenu','Position',[150 hPosition1 120 20],...
    'String',{'NorthEast', 'SouthEast', 'SouthWest', 'NorthWest'},'Value',3,...
    'Tag','popupmenu_imageScaleBarLocation','Enable',scaleBarStatus,...
    'Callback',@(h,event) setScaleBar(guidata(h),'imageScaleBar'));

hPosition1=hPosition1+30;
uicontrol(imagePanel,'Style','text','Position',[20 hPosition1-2 50 20],...
    'String','Cmin','HorizontalAlignment','left');
uicontrol(imagePanel,'Style','edit','Position',[70 hPosition1 50 20],...
    'String','','BackgroundColor','white','Tag','edit_cmin',...
    'Callback',@(h,event) setClim(guidata(h)));
uicontrol(imagePanel,'Style','text','Position',[150 hPosition1-2 50 20],...
    'String','Cmax','HorizontalAlignment','left');
uicontrol(imagePanel,'Style','edit','Position',[200 hPosition1 50 20],...
    'String','','BackgroundColor','white','Tag','edit_cmax',...
    'Callback',@(h,event) setClim(guidata(h)));

hPosition1=hPosition1+20;
uicontrol(imagePanel,'Style','checkbox',...
    'Position',[10 hPosition1 200 20],'Tag','checkbox_autoscale',...
    'String','Autoscale','HorizontalAlignment','left','FontWeight','bold',...
    'Callback',@(h,event) setClim(guidata(h)));

hPosition1=hPosition1+50;
% Create image processes controls
nProc = numel(imageProc);
for iProc=nProc:-1:1;
    output=imageProc{iProc}.getDrawableOutput;
    validChan = imageProc{iProc}.checkChannelOutput;
    validOutput = find(strcmp({output.type},'image'));
    for iOutput=validOutput(end:-1:1)
        createOutputText(imagePanel,imageProcId(iProc),iOutput,hPosition1,output(iOutput).name);
        arrayfun(@(x) createChannelButton(imagePanel,imageProcId(iProc),iOutput,x,hPosition1),...
            find(validChan));
        hPosition1=hPosition1+20;
    end
    createProcText(imagePanel,imageProcId(iProc),iOutput,hPosition1,imageProc{iProc}.getName);
    hPosition1=hPosition1+20;
end

% Create channels controls
hPosition1=hPosition1+10;
uicontrol(imagePanel,'Style','radio','Position',[10 hPosition1 200 20],...
    'Tag','radiobutton_channels','String',' Channels','Value',1,...
    'HorizontalAlignment','left','FontWeight','bold');
arrayfun(@(i) uicontrol(imagePanel,'Style','checkbox',...
    'Position',[200+30*i hPosition1 20 20],...
    'Tag',['checkbox_channel' num2str(i)],'Value',i<3,...
    'Callback',@(h,event) redrawChannel(h,guidata(h))),...
    1:numel(userData.MD.channels_));
hPosition1=hPosition1+20;
arrayfun(@(i) uicontrol(imagePanel,'Style','text',...
    'Position',[200+30*i hPosition1 20 20],...
    'Tag',['text_channel' num2str(i)],'String',i),...
    1:numel(userData.MD.channels_));

%% Create overlay panel
overlayPanel = uipanel(mainFig,'Position',[1/2 0 1/2 1],...
    'Title','Overlay','BackgroundColor',get(0,'defaultUicontrolBackgroundColor'),...
    'Units','pixels','Tag','uipanel_overlay');

% Create image option controls
hPosition2=10;
uicontrol(overlayPanel,'Style','checkbox',...
    'Position',[20 hPosition2 100 20],'Tag','checkbox_vectorFieldScaleBar',...
    'String',' Scale bar','HorizontalAlignment','left',...
    'Enable',scaleBarStatus,...
    'Callback',@(h,event) setScaleBar(guidata(h),'vectorFieldScaleBar'));
uicontrol(overlayPanel,'Style','edit','Position',[120 hPosition2 50 20],...
    'String','1000','BackgroundColor','white','Tag','edit_vectorFieldScaleBar',...
    'Enable',scaleBarStatus,...
    'Callback',@(h,event) setScaleBar(guidata(h),'vectorFieldScaleBar'));
uicontrol(overlayPanel,'Style','text','Position',[180 hPosition2 70 20],...
    'String','nm/min','HorizontalAlignment','left');

hPosition2=hPosition2+30;
uicontrol(overlayPanel,'Style','text',...
    'Position',[20 hPosition2 100 20],'Tag','text_vectorFieldScale',...
    'String',' Vector scale','HorizontalAlignment','left');
uicontrol(overlayPanel,'Style','edit','Position',[120 hPosition2 50 20],...
    'String','10','BackgroundColor','white','Tag','edit_vectorFieldScale',...
    'Callback',@(h,event) redrawOverlays(guidata(h)));

hPosition2=hPosition2+20;
uicontrol(overlayPanel,'Style','text',...
    'Position',[10 hPosition2 200 20],'Tag','text_vectorFieldOptions',...
    'String','Vector field options','HorizontalAlignment','left','FontWeight','bold');
hPosition2=hPosition2+50;


nProc = numel(overlayProc);
for iProc=nProc:-1:1;
    output=overlayProc{iProc}.getDrawableOutput;
    validChan = overlayProc{iProc}.checkChannelOutput;
    validOutput = find(strcmp({output.type},'overlay'));
    for iOutput=validOutput(end:-1:1)
        createOutputText(overlayPanel,overlayProcId(iProc),iOutput,hPosition2,output(iOutput).name);
        arrayfun(@(x) createChannelBox(overlayPanel,overlayProcId(iProc),iOutput,x,hPosition2),...
            find(validChan));
        hPosition2=hPosition2+20;
    end
    createProcText(overlayPanel,overlayProcId(iProc),iOutput,hPosition2,overlayProc{iProc}.getName);
    hPosition2=hPosition2+20;
end
arrayfun(@(i) uicontrol(overlayPanel,'Style','text',...
    'Position',[200+30*i hPosition2 20 20],...
    'Tag',['text_channel' num2str(i)],'String',i),...
    1:numel(userData.MD.channels_));

%% Resize panels
% Get image panel size
a=get(get(imagePanel,'Children'),'Position');
P=vertcat(a{:});
imagePanelSize = [max(P(:,1)+P(:,3))+20 max(P(:,2)+P(:,4))+20];
% Get overlay panel size
a=get(get(overlayPanel,'Children'),'Position');
P=vertcat(a{:});
overlayPanelSize = [max(P(:,1)+P(:,3))+20 max(P(:,2)+P(:,4))+20];

panelsLength = imagePanelSize(1)+overlayPanelSize(1)+10;
panelsHeight = max(imagePanelSize(2),overlayPanelSize(2));

% Resize panel
set(imagePanel,'Position',[10 panelsHeight-imagePanelSize(2)+10 ...
    imagePanelSize(1) imagePanelSize(2)],...
    'SelectionChangeFcn',@(h,event) redrawImage(guidata(h)))
set(overlayPanel,'Position',[10+imagePanelSize(1)+10 panelsHeight-overlayPanelSize(2)+10 ...
    overlayPanelSize(1) overlayPanelSize(2)])

%% Create movie panel

parentPanel = uipanel(mainFig,...
    'Title','','BackgroundColor',get(0,'defaultUicontrolBackgroundColor'),...
    'Units','pixels','Tag','uipanel_movie','BorderType','none');
set(parentPanel,'Position',[10 panelsHeight+10 panelsLength 100]);

uicontrol(parentPanel,'Style','text','Position',[10 10 50 15],...
    'String','Frame','Tag','text_frame','HorizontalAlignment','left');
uicontrol(parentPanel,'Style','edit','Position',[70 10 30 20],...
    'String','1','Tag','edit_frame','BackgroundColor','white',...
    'HorizontalAlignment','left',...
    'Callback',@(h,event) redrawScene(h,guidata(h)));
uicontrol(parentPanel,'Style','text','Position',[100 10 40 15],...
    'HorizontalAlignment','left',...
    'String',['/' num2str(userData.MD.nFrames_)],'Tag','text_frameMax');

uicontrol(parentPanel,'Style','slider',...
    'Position',[150 10 panelsLength-160 20],...
    'Value',1,'Min',1,'Max',userData.MD.nFrames_,...
    'SliderStep',[1/double(userData.MD.nFrames_)  5/double(userData.MD.nFrames_)],...
    'Tag','slider_frame','BackgroundColor','white',...
    'Callback',@(h,event) redrawScene(h,guidata(h)));

uicontrol(parentPanel,'Style','text','Position',[10 40 40 20],...
    'String','Movie','Tag','text_movie');
uicontrol(parentPanel,'Style','edit','Position',[60 40 panelsLength-70 20],...
    'String',[userData.MD.movieDataPath_ filesep userData.MD.movieDataFileName_],...
    'HorizontalAlignment','left','BackgroundColor','white','Tag','edit_movie');

uicontrol(parentPanel,'Style','text','Position',[10 70 panelsLength 20],...
    'String',userfcn_softwareConfig(),'Tag','text_copyright',...
    'HorizontalAlignment','left');


%% Resize panels and figure
sz=get(0,'ScreenSize');
figWidth = panelsLength+20;
figHeight = panelsHeight+110;
set(mainFig,'Position',[sz(3)/50 (sz(4)-figHeight)/2 figWidth figHeight]);

% Auto select input process
if ismember(ip.Results.procId,validProcId)
    for i=ip.Results.procId
        h=findobj(mainFig,'-regexp','Tag',['(\w)_process' ...
            num2str(i)  '_output1_channel(\d+)']);
        set(h,'Value',1);
    end
end

% Update handles structure
handles = guihandles(mainFig);
guidata(handles.figure1, handles);

userData.drawFig=-1;
set(handles.figure1,'UserData',userData);

% Update the image and overlays
redrawScene(handles.figure1, handles);
% playMovie(handles.figure1,handles)

function playMovie(hObject,handles)

userData = get(handles.figure1, 'UserData');
nf = userData.MD.nFrames_;
nf=10;
fmt = ['%0' num2str(ceil(log10(nf))) 'd'];
fpath = [userData.MD.movieDataPath_ filesep 'frames' filesep];
mkClrDir(fpath);
mpath = [userData.MD.movieDataPath_ filesep 'movie' filesep];
mkClrDir(mpath);
fprintf('Generating movie frames:     ');
for f=1:nf
    set(handles.slider_frame, 'Value',f);
    redrawScene(hObject, handles);
    drawnow;
    print(userData.drawFig, '-dpng', '-loose', ['-r' num2str(1*72)], [fpath 'frame' num2str(f, fmt) '.png']);
    fprintf('\b\b\b\b%3d%%', round(100*f/(nf)));
end
fprintf('\n');

% Generate movie
fprintf('Generating movie... ');
fr = num2str(15);
cmd = ['ffmpeg -loglevel quiet -y -r ' fr ' -i ' fpath 'frame' fmt '.png' ' -r ' fr ' -b 50000k -bt 20000k ' mpath 'movie.mp4 > /dev/null 2>&1' ];
system(cmd);
fprintf('done.\n');

function redrawScene(hObject, handles)

userData = get(handles.figure1, 'UserData');
% Retrieve the value of the selected image
if strcmp(get(hObject,'Tag'),'edit_frame')
    frameNumber = str2double(get(handles.edit_frame, 'String'));
else
    frameNumber = get(handles.slider_frame, 'Value');
end
frameNumber=round(frameNumber);
frameNumber = min(max(frameNumber,1),userData.MD.nFrames_);

% Set the slider and editboxes values
set(handles.edit_frame,'String',frameNumber);
set(handles.slider_frame,'Value',frameNumber);

% Update the image and overlays
redrawImage(handles);
redrawOverlays(handles);

function createFigure(handles)
userData = get(handles.figure1,'UserData');
%Create a figure
sz=get(0,'ScreenSize');
ratios = [sz(3)/userData.MD.imSize_(2) sz(4)/userData.MD.imSize_(1)];
nx=.6*min(ratios)*userData.MD.imSize_(2);
ny=.6*min(ratios)*userData.MD.imSize_(1);
userData.drawFig = figure('Position',[sz(3)*.2 sz(4)*.2 nx ny],...
    'Name','Figure','NumberTitle','off');

% figure options for movie export
% iptsetpref('ImshowBorder','tight');
set(userData.drawFig, 'InvertHardcopy', 'off');
set(userData.drawFig, 'PaperUnits', 'Points');
set(userData.drawFig, 'PaperSize', [nx ny]);
set(userData.drawFig, 'PaperPosition', [0 0 nx ny]); % very important
%  set(userData.drawFig,'DefaultLineLineSmoothing','on');
% set(userData.drawFig,'DefaultPatchLineSmoothing','on');

%Create the associate axes
axes('Parent',userData.drawFig,'XLim',[0 userData.MD.imSize_(2)],...
    'YLim',[0 userData.MD.imSize_(1)],'Position',[0.05 0.05 .9 .9]);
set(handles.figure1,'UserData',userData);


function redrawChannel(hObject,handles)
% Callback for channels checkboxes to avoid 0 or more than 4 channels
channelBoxes = findobj(handles.figure1,'-regexp','Tag','checkbox_channel*');
chanList=find(arrayfun(@(x)get(x,'Value'),channelBoxes));
if numel(chanList)==0
    set(hObject,'Value',1);
elseif numel(chanList)>4
   set(hObject,'Value',0); 
end

redrawImage(handles)

function setScaleBar(handles,type)
% Remove existing scalebar of given type if any
h=findobj('Tag',type);
if ~isempty(h), delete(h); end

% If checked, adds a new scalebar using the width as a label input
userData=get(handles.figure1,'UserData');
if ~get(handles.(['checkbox_' type]),'Value') || ~ishandle(userData.drawFig),
    return 
end
figure(userData.drawFig)
scale = str2double(get(handles.(['edit_' type]),'String'));
if strcmp(type,'imageScaleBar')
    width = scale *1000/userData.MD.pixelSize_;
    if get(handles.(['checkbox_' type 'Label']),'Value')
        label = [num2str(scale) ' \mum'];
    else
        label='';
    end
    props=get(handles.(['popupmenu_' type 'Location']),{'String','Value'});
    location=props{1}{props{2}};
else
    displayScale = str2double(get(handles.edit_vectorFieldScale,'String'));
    width = scale*displayScale/(userData.MD.pixelSize_/userData.MD.timeInterval_*60);
    label= [num2str(scale) ' nm/min'];
    location='SouthEast';
end
hScaleBar = plotScaleBar(width,'Label',label,'Location',location);
set(hScaleBar,'Tag',type);


function setTimeStamp(handles)
% Remove existing scalebar of given type if any
h=findobj('Tag','timeStamp');
if ~isempty(h), delete(h); end

% If checked, adds a new scalebar using the width as a label input
userData=get(handles.figure1,'UserData');
if ~get(handles.checkbox_timeStamp,'Value') || ~ishandle(userData.drawFig),
    return 
end
figure(userData.drawFig)
frameNr=get(handles.slider_frame,'Value');
width = userData.MD.imSize_(2)/20;
time= (frameNr-1)*userData.MD.timeInterval_;
p=sec2struct(time);
props=get(handles.popupmenu_timeStampLocation,{'String','Value'});
location=props{1}{props{2}};
hTimeStamp = plotScaleBar(width,'Label',p.str,'Location',location);
set(hTimeStamp,'Tag','timeStamp');
delete(hTimeStamp(1))

function setClim(handles)
userData=get(handles.figure1,'UserData');
imageTag = get(get(handles.uipanel_image,'SelectedObject'),'Tag');

clim=[str2double(get(handles.edit_cmin,'String')) ...
    str2double(get(handles.edit_cmax,'String'))];

if strcmp(imageTag,'radiobutton_channels')
    channelBoxes = findobj(handles.figure1,'-regexp','Tag','checkbox_channel*');
    chanList=find(arrayfun(@(x)get(x,'Value'),channelBoxes));
    userData.MD.channels_(chanList(1)).displayMethod_.CLim=cLim;
else
    % Retrieve the id, process nr and channel nr of the selected imageProc
    tokens = regexp(imageTag,'radiobutton_process(\d+)_output(\d+)_channel(\d+)','tokens');
    procId=str2double(tokens{1}{1});
    iOutput = str2double(tokens{1}{2});
    iChan = str2double(tokens{1}{3});
    userData.MD.processes_{procId}.displayMethod_{iOutput,iChan}.CLim=clim;
end

if ishandle(userData.drawFig)
    child=get(userData.drawFig,'Children');
    set(child(strcmp(get(child,'Type'),'axes')),'Clim',clim);
end


function redrawImage(handles)
userData=get(handles.figure1,'UserData');
frameNr=get(handles.slider_frame,'Value');

imageTag = get(get(handles.uipanel_image,'SelectedObject'),'Tag');

if ~ishandle(userData.drawFig), 
    createFigure(handles);
else
    figure(userData.drawFig);
end
channelBoxes = findobj(handles.figure1,'-regexp','Tag','checkbox_channel*');
if strcmp(imageTag,'radiobutton_channels')
    set(channelBoxes,'Enable','on');
    chanList=find(arrayfun(@(x)get(x,'Value'),channelBoxes));
    userData.MD.channels_(chanList).draw(frameNr);
    clim=userData.MD.channels_(chanList(1)).displayMethod_.CLim;
else
    set(channelBoxes,'Enable','off');
    % Retrieve the id, process nr and channel nr of the selected imageProc
    tokens = regexp(imageTag,'radiobutton_process(\d+)_output(\d+)_channel(\d+)','tokens');
    procId=str2double(tokens{1}{1});
    outputList = userData.MD.processes_{procId}.getDrawableOutput;
    iOutput = str2double(tokens{1}{2});
    output = outputList(iOutput).var;
    iChan = str2double(tokens{1}{3});
    userData.MD.processes_{procId}.draw(iChan,frameNr,'output',output);
    clim=userData.MD.processes_{procId}.displayMethod_{iOutput,iChan}.CLim;
end

% Set the autoscale
set(handles.checkbox_autoscale,'Value',isempty(clim));
if isempty(clim)
    set([handles.edit_cmin handles.edit_cmax],'Enable','off');
else
    set(handles.edit_cmin,'Enable','on','String',clim(1));
    set(handles.edit_cmax,'Enable','on','String',clim(2));
end

% Reset the scaleBar
setScaleBar(handles,'imageScaleBar'); 
setTimeStamp(handles); 

function redrawOverlays(handles)
overlayBoxes = findobj(handles.uipanel_overlay,'-regexp','Tag','checkbox_process*');
checkedBoxes = logical(arrayfun(@(x) get(x,'Value'),overlayBoxes));
overlayTags=arrayfun(@(x) get(x,'Tag'),overlayBoxes(checkedBoxes),...
    'UniformOutput',false);
for i=1:numel(overlayTags)
    redrawOverlay(handles.(overlayTags{i}),handles)
end

% Reset the scaleBar
if get(handles.checkbox_vectorFieldScaleBar,'Value'), 
    setScaleBar(handles,'vectorFieldScaleBar');
end
    
function redrawOverlay(hObject,handles)
userData=get(handles.figure1,'UserData');
frameNr=get(handles.slider_frame,'Value');

overlayTag = get(hObject,'Tag');
if ishandle(userData.drawFig), 
    figure(userData.drawFig); 
else
    redrawScene(hObject, handles); return;
end
 % Retrieve the id, process nr and channel nr of the selected imageProc
tokens = regexp(overlayTag,'checkbox_process(\d+)_output(\d+)_channel(\d+)','tokens');
procId=str2double(tokens{1}{1});
outputList = userData.MD.processes_{procId}.getDrawableOutput;
iOutput = str2double(tokens{1}{2});
output = outputList(iOutput).var;
iChan = str2double(tokens{1}{3});
if get(hObject,'Value')
    userData.MD.processes_{procId}.draw(iChan,frameNr,'output',output,...
        'scale',str2double(get(handles.edit_vectorFieldScale,'String')));
else
    h=findobj('Tag',[userData.MD.processes_{procId}.getName '_channel'...
        num2str(iChan) '_output' num2str(iOutput)]);
    if ~isempty(h), delete(h); end
end

function deleteViewer(handles)

userData=get(handles.figure1,'UserData');
if isfield(userData, 'drawFig') && ishandle(userData.drawFig), 
    delete(userData.drawFig); 
end