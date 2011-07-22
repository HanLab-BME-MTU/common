function mainFig = movieDataViewer(varargin)

ip = inputParser;
ip.addRequired('MD',@(x) isa(x,'MovieData'));
ip.addOptional('procId',0,@isscalar);
ip.parse(varargin{:});

% Chek
h=findobj(0,'Name','MovieDataViewer');
if ~isempty(h), delete(h); end
mainFig=figure('Name','MovieDataViewer','Position',[0 0 200 200],...
    'NumberTitle','off','Tag','figure1','Toolbar','none','MenuBar','none',...
    'Color',get(0,'defaultUicontrolBackgroundColor'),...
    'DeleteFcn', @(h,event) deleteViewer(h,event,guidata(h)));
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
    'Callback',@(h,event) checkOverlay(h,event,guidata(h)));

%% Create image panel
imagePanel = uibuttongroup(mainFig,'Position',[0 0 1/2 1],...
    'Title','Image','BackgroundColor',get(0,'defaultUicontrolBackgroundColor'),...
    'Units','pixels','Tag','uipanel_image');

% Create image processes controls
nProc = numel(imageProc);
hPosition1=10;
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
    'Callback',@(h,event) checkChannel(h,event,guidata(h))),...
    1:numel(userData.MD.channels_));
hPosition1=hPosition1+20;
arrayfun(@(i) uicontrol(imagePanel,'Style','text',...
    'Position',[200+30*i hPosition1 20 20],...
    'Tag',['text_channel' num2str(i)],'String',i),...
    1:numel(userData.MD.channels_));
a=get(get(imagePanel,'Children'),'Position');
P=vertcat(a{:});
imagePanelSize = [max(P(:,1)+P(:,3))+20 max(P(:,2)+P(:,4))+20];
% Resize panel
set(imagePanel,'Position',[10 10 imagePanelSize(1) imagePanelSize(2)],...
    'SelectionChangeFcn',@(h,event) redrawImage(h,event,guidata(h)))

%% Create overlay panel
overlayPanel = uipanel(mainFig,'Position',[1/2 0 1/2 1],...
    'Title','Overlay','BackgroundColor',get(0,'defaultUicontrolBackgroundColor'),...
    'Units','pixels','Tag','uipanel_overlay');
hPosition2=10;
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
% Resize panel
set(imagePanel,'Position',[10 10 imagePanelSize(1) imagePanelSize(2)],...
    'SelectionChangeFcn',@(h,event) redrawImage(h,event,guidata(h)))
% Get overlay panel size
a=get(get(overlayPanel,'Children'),'Position');
P=vertcat(a{:});
overlayPanelSize = [max(P(:,1)+P(:,3))+20 max(P(:,2)+P(:,4))+20];

panelsLength = imagePanelSize(1)+overlayPanelSize(1)+10;
panelsHeight = max(imagePanelSize(2),overlayPanelSize(2));

% Resize panel
set(imagePanel,'Position',[10 panelsHeight-imagePanelSize(2)+10 ...
    imagePanelSize(1) imagePanelSize(2)],...
    'SelectionChangeFcn',@(h,event) redrawImage(h,event,guidata(h)))
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
    'Callback',@(h,event) frameEdition(h,event,guidata(h)));
uicontrol(parentPanel,'Style','text','Position',[100 10 40 15],...
    'HorizontalAlignment','left',...
    'String',['/' num2str(userData.MD.nFrames_)],'Tag','text_frameMax');

uicontrol(parentPanel,'Style','slider',...
    'Position',[150 10 panelsLength-160 20],...
    'Value',1,'Min',1,'Max',userData.MD.nFrames_,...
    'SliderStep',[1/double(userData.MD.nFrames_)  5/double(userData.MD.nFrames_)],...
    'Tag','slider_frame','BackgroundColor','white',...
    'Callback',@(h,event) frameEdition(h,event,guidata(h)));

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
set(mainFig,'Position',[100 sz(4)-100   panelsLength+20 panelsHeight+100+10]);

% Auto select input process
if ismember(ip.Results.procId,validProcId)
    h=findobj(mainFig,'-regexp','Tag',['(\w)_process' ...
        num2str(ip.Results.procId)  '_output1_channel(\d+)']);
    set(h,'Value',1);
end

% Update handles structure
handles = guihandles(mainFig);
guidata(handles.figure1, handles);

userData.drawFig=-1;
set(handles.figure1,'UserData',userData);

% Update the image and overlays
redrawImage(handles.figure1, [], handles);
redrawOverlays(handles.figure1, [], handles);

% --- Executes on slider movement.
function frameEdition(hObject, eventdata, handles)

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

% Update the image ad overlays
redrawImage(hObject, eventdata, handles);
redrawOverlays(hObject, eventdata, handles);


function checkChannel(hObject,event,handles)
% Callback of channels checkboxes to avoid 0 or more than 4 channels
channelBoxes = findobj(handles.figure1,'-regexp','Tag','checkbox_channel*');
chanList=find(arrayfun(@(x)get(x,'Value'),channelBoxes));
if numel(chanList)==0
    set(hObject,'Value',1);
elseif numel(chanList)>4
   set(hObject,'Value',0); 
end

redrawImage(hObject,event,handles)

function redrawImage(hObject, eventdata, handles)
userData=get(handles.figure1,'UserData');
frameNr=get(handles.slider_frame,'Value');

imageTag = get(get(handles.uipanel_image,'SelectedObject'),'Tag');

if ishandle(userData.drawFig), 
    figure(userData.drawFig); 
else 
    %Create a figure
    sz=get(0,'ScreenSize');
    ratios = [sz(3)/userData.MD.imSize_(2) sz(4)/userData.MD.imSize_(1)];
    userData.drawFig = figure('Position',[sz(3)*.2 sz(4)*.2 ...
        .6*min(ratios)*userData.MD.imSize_(2) .6*min(ratios)*userData.MD.imSize_(1)]);
    
    %Create the associate axes
    axes('Parent',userData.drawFig,'XLim',[0 userData.MD.imSize_(2)],'YLim',[0 userData.MD.imSize_(1)],...
        'Position',[0.05 0.05 .9 .9]);
end
set(handles.figure1,'UserData',userData);

channelBoxes = findobj(handles.figure1,'-regexp','Tag','checkbox_channel*');
if strcmp(imageTag,'radiobutton_channels')
    set(channelBoxes,'Enable','on');
    chanList=logical(arrayfun(@(x)get(x,'Value'),channelBoxes));
    userData.MD.channels_(chanList).draw(frameNr);
else
    set(channelBoxes,'Enable','off');
    % Retrieve the id, process nr and channel nr of the selected imageProc
    tokens = regexp(imageTag,'radiobutton_process(\d+)_output(\d+)_channel(\d+)','tokens');
    procId=str2double(tokens{1}{1});
    outputList = userData.MD.processes_{procId}.getDrawableOutput;
    output = outputList(str2double(tokens{1}{2})).var;
    iChan = str2double(tokens{1}{3});
    userData.MD.processes_{procId}.draw(iChan,frameNr,'output',output);
end

function checkOverlay(hObject, eventdata, handles)
userData=get(handles.figure1,'UserData');
frameNr=get(handles.slider_frame,'Value');

overlayTag = get(hObject,'Tag');

if ishandle(userData.drawFig), figure(userData.drawFig); end
 % Retrieve the id, process nr and channel nr of the selected imageProc
tokens = regexp(overlayTag,'checkbox_process(\d+)_output(\d+)_channel(\d+)','tokens');
procId=str2double(tokens{1}{1});
outputList = userData.MD.processes_{procId}.getDrawableOutput;
iOutput = str2double(tokens{1}{2});
output = outputList(iOutput).var;
iChan = str2double(tokens{1}{3});
if get(hObject,'Value')
    userData.MD.processes_{procId}.draw(iChan,frameNr,'output',output);
else
    h=findobj('Tag',[userData.MD.processes_{procId}.getName '_channel'...
        num2str(iChan) '_output' num2str(iOutput)]);
    if ~isempty(h), delete(h); end
end

function redrawOverlays(hObject, eventdata, handles)
userData=get(handles.figure1,'UserData');
frameNr=get(handles.slider_frame,'Value');

figure(userData.drawFig);
overlayBoxes = findobj(handles.uipanel_overlay,'Style','checkbox');
checkedBoxes = logical(arrayfun(@(x) get(x,'Value'),overlayBoxes));
overlayTags=arrayfun(@(x) get(x,'Tag'),overlayBoxes(checkedBoxes),...
    'UniformOutput',false);
for i=1:numel(overlayTags)
    overlayTag=overlayTags{i};
    tokens = regexp(overlayTag,'checkbox_process(\d+)_output(\d+)_channel(\d+)','tokens');
    procId=str2double(tokens{1}{1});
    outputList = userData.MD.processes_{procId}.getDrawableOutput;
    output = outputList(str2double(tokens{1}{2})).var;
    iChan = str2double(tokens{1}{3});

    userData.MD.processes_{procId}.draw(iChan,frameNr,'output',output);
end


% --- Executes during object deletion, before destroying properties.
function deleteViewer(hObject, eventdata, handles)

userData=get(handles.figure1,'UserData');
if isfield(userData, 'drawFig') && ishandle(userData.drawFig), 
    delete(userData.drawFig); 
end