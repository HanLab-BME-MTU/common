function mainFig = movieViewer(MO,varargin)

ip = inputParser;
ip.addRequired('MO',@(x) isa(x,'MovieObject'));
ip.addOptional('procId',[],@isnumeric);
ip.addParamValue('movieIndex',0,@isscalar);
ip.parse(MO,varargin{:});

% Chek
h=findobj(0,'Name','Viewer');
if ~isempty(h), delete(h); end
mainFig=figure('Name','Viewer','Position',[0 0 200 200],...
    'NumberTitle','off','Tag','figure1','Toolbar','none','MenuBar','none',...
    'Color',get(0,'defaultUicontrolBackgroundColor'),'Resize','off',...
    'DeleteFcn', @(h,event) deleteViewer());
userData=get(mainFig,'UserData');

if isa(ip.Results.MO,'MovieList')
    userData.ML=ip.Results.MO;
    userData.movieIndex=ip.Results.movieIndex;
    if userData.movieIndex~=0
        userData.MO=ip.Results.MO.movies_{userData.movieIndex};
    else
         userData.MO=ip.Results.MO;
    end
        
%     userData.MO=MO.movies_{userData.movieIndex};
    userData.procId = ip.Results.procId;
    if ~isempty(ip.Results.procId)
        procId = userData.MO.getProcessIndex(class(userData.ML.processes_{ip.Results.procId}));
    else
        procId = ip.Results.procId;
    end
else
    userData.MO=ip.Results.MO;
%     userData.MO=ip.Results.MO;
    procId=ip.Results.procId;
end

% Classify movieData processes by type (image, overlay, movie overlay or
% graph)
validProcId= find(cellfun(@(x) ismember('getDrawableOutput',methods(x)) &...
    x.success_,userData.MO.processes_));
validProc=userData.MO.processes_(validProcId);

getOutputType = @(type) cellfun(@(x) any(~cellfun(@isempty,regexp({x.getDrawableOutput.type},type,'once','start'))),...
    validProc);

isImageProc =getOutputType('image');
imageProc=validProc(isImageProc);
imageProcId = validProcId(isImageProc);
isOverlayProc =getOutputType('[oO]verlay');
overlayProc=validProc(isOverlayProc);
overlayProcId = validProcId(isOverlayProc);
isGraphProc =getOutputType('[gG]raph');
graphProc=validProc(isGraphProc);
graphProcId = validProcId(isGraphProc);

% Create series of anonymous function to generate process controls
createProcText= @(panel,i,j,pos,name) uicontrol(panel,'Style','text',...
    'Position',[10 pos 250 20],'Tag',['text_process' num2str(i)],...
    'String',name,'HorizontalAlignment','left','FontWeight','bold');
createOutputText= @(panel,i,j,pos,text) uicontrol(panel,'Style','text',...
    'Position',[40 pos 200 20],'Tag',['text_process' num2str(i) '_output'...
    num2str(j)],'String',text,'HorizontalAlignment','left');
createProcButton= @(panel,i,j,k,pos) uicontrol(panel,'Style','radio',...
    'Position',[200+30*k pos 20 20],'Tag',['radiobutton_process' num2str(i) '_output'...
    num2str(j) '_channel' num2str(k)]);
createChannelBox= @(panel,i,j,k,pos,varargin) uicontrol(panel,'Style','checkbox',...
    'Position',[200+30*k pos 20 20],'Tag',['checkbox_process' num2str(i) '_output'...
    num2str(j) '_channel' num2str(k)],varargin{:});
createMovieBox= @(panel,i,j,pos,name,varargin) uicontrol(panel,'Style','checkbox',...
    'Position',[40 pos 200 25],'Tag',['checkbox_process' num2str(i) '_output'...
    num2str(j)],'String',[' ' name],varargin{:});
createInputBox= @(panel,i,j,k,pos,name,varargin) uicontrol(panel,'Style','checkbox',...
    'Position',[40 pos 200 25],'Tag',['checkbox_process' num2str(i) '_output'...
    num2str(j) '_input' num2str(k)],'String',[' ' name],varargin{:});
createInputInputBox= @(panel,i,j,k,l,pos,varargin) uicontrol(panel,'Style','checkbox',...
    'Position',[200+30*l pos 20 20],'Tag',['checkbox_process' num2str(i) '_output'...
    num2str(j) '_input' num2str(k) '_input' num2str(l)],varargin{:});


%% Image panel creation
if isa(userData.MO,'MovieData')
    imagePanel = uibuttongroup(mainFig,'Position',[0 0 1/2 1],...
        'Title','Image','BackgroundColor',get(0,'defaultUicontrolBackgroundColor'),...
        'Units','pixels','Tag','uipanel_image');
    
    % First create image option (timestamp, scalebar, image scaling)
    % Timestamp
    hPosition1=10;
    if isempty(userData.MO.timeInterval_),
        timeStampStatus = 'off';
    else
        timeStampStatus = 'on';
    end
    uicontrol(imagePanel,'Style','checkbox',...
        'Position',[10 hPosition1 200 20],'Tag','checkbox_timeStamp',...
        'String',' Time stamp','HorizontalAlignment','left',...
        'Enable',timeStampStatus,'Callback',@(h,event) setTimeStamp(guidata(h)));
    uicontrol(imagePanel,'Style','popupmenu','Position',[130 hPosition1 120 20],...
        'String',{'NorthEast', 'SouthEast', 'SouthWest', 'NorthWest'},'Value',4,...
        'Tag','popupmenu_timeStampLocation','Enable',timeStampStatus,...
        'Callback',@(h,event) setTimeStamp(guidata(h)));
    
    % Scalebar
    hPosition1=hPosition1+30;
    if isempty(userData.MO.pixelSize_),
        scaleBarStatus = 'off';
    else
        scaleBarStatus = 'on';
    end
    uicontrol(imagePanel,'Style','edit','Position',[30 hPosition1 50 20],...
        'String','1','BackgroundColor','white','Tag','edit_imageScaleBar',...
        'Enable',scaleBarStatus,...
        'Callback',@(h,event) setScaleBar(guidata(h),'imageScaleBar'));
    uicontrol(imagePanel,'Style','text','Position',[85 hPosition1-2 70 20],...
        'String','microns','HorizontalAlignment','left');
    uicontrol(imagePanel,'Style','checkbox',...
        'Position',[150 hPosition1 100 20],'Tag','checkbox_imageScaleBarLabel',...
        'String',' Show label','HorizontalAlignment','left',...
        'Enable',scaleBarStatus,...
        'Callback',@(h,event) setScaleBar(guidata(h),'imageScaleBar'));
    
    hPosition1=hPosition1+30;
    uicontrol(imagePanel,'Style','checkbox',...
        'Position',[10 hPosition1 200 20],'Tag','checkbox_imageScaleBar',...
        'String',' Scalebar','HorizontalAlignment','left',...
        'Enable',scaleBarStatus,...
        'Callback',@(h,event) setScaleBar(guidata(h),'imageScaleBar'));
    uicontrol(imagePanel,'Style','popupmenu','Position',[130 hPosition1 120 20],...
        'String',{'NorthEast', 'SouthEast', 'SouthWest', 'NorthWest'},'Value',3,...
        'Tag','popupmenu_imageScaleBarLocation','Enable',scaleBarStatus,...
        'Callback',@(h,event) setScaleBar(guidata(h),'imageScaleBar'));
    
    % Colormap control
    hPosition1=hPosition1+30;
    uicontrol(imagePanel,'Style','text','Position',[20 hPosition1-2 100 20],...
        'String','Color limits','HorizontalAlignment','left');
    uicontrol(imagePanel,'Style','edit','Position',[150 hPosition1 50 20],...
        'String','','BackgroundColor','white','Tag','edit_cmin',...
        'Callback',@(h,event) setCLim(guidata(h)));
    uicontrol(imagePanel,'Style','edit','Position',[200 hPosition1 50 20],...
        'String','','BackgroundColor','white','Tag','edit_cmax',...
        'Callback',@(h,event) setCLim(guidata(h)));
    
    hPosition1=hPosition1+30;
    uicontrol(imagePanel,'Style','checkbox',...
        'Position',[10 hPosition1 120 20],'Tag','checkbox_colorbar',...
        'String',' Colorbar','HorizontalAlignment','left',...
        'Callback',@(h,event) setColorbar(guidata(h)));
    
    uicontrol(imagePanel,'Style','text','Position',[120 hPosition1-2 80 20],...
        'String','Colormap','HorizontalAlignment','left');
    uicontrol(imagePanel,'Style','popupmenu',...
        'Position',[200 hPosition1 80 20],'Tag','popupmenu_colormap',...
        'String',{'Gray','Jet','HSV'},'Value',1,...
        'HorizontalAlignment','left','Callback',@(h,event) setColormap(guidata(h)));
    
    hPosition1=hPosition1+20;
    uicontrol(imagePanel,'Style','text','Position',[10 hPosition1 200 20],...
        'String','Image options','HorizontalAlignment','left','FontWeight','bold');
    
    
    
    % Create controls for switching between process image output
    hPosition1=hPosition1+50;
    nProc = numel(imageProc);
    for iProc=nProc:-1:1;
        output=imageProc{iProc}.getDrawableOutput;
        validChan = imageProc{iProc}.checkChannelOutput;
        validOutput = find(strcmp({output.type},'image'));
        for iOutput=validOutput(end:-1:1)
            createOutputText(imagePanel,imageProcId(iProc),iOutput,hPosition1,output(iOutput).name);
            arrayfun(@(x) createProcButton(imagePanel,imageProcId(iProc),iOutput,x,hPosition1),...
                find(validChan));
            hPosition1=hPosition1+20;
        end
        createProcText(imagePanel,imageProcId(iProc),iOutput,hPosition1,imageProc{iProc}.getName);
        hPosition1=hPosition1+20;
    end
    
    % Create controls for selecting channels (raw image)
    hPosition1=hPosition1+10;
    uicontrol(imagePanel,'Style','radio','Position',[10 hPosition1 200 20],...
        'Tag','radiobutton_channels','String',' Raw image','Value',1,...
        'HorizontalAlignment','left','FontWeight','bold');
    arrayfun(@(i) uicontrol(imagePanel,'Style','checkbox',...
        'Position',[200+30*i hPosition1 20 20],...
        'Tag',['checkbox_channel' num2str(i)],'Value',i<4,...
        'Callback',@(h,event) redrawChannel(h,guidata(h))),...
        1:numel(userData.MO.channels_));
    
    hPosition1=hPosition1+20;
    uicontrol(imagePanel,'Style','text','Position',[120 hPosition1 100 20],...
        'Tag','text_channels','String','Channels');
    arrayfun(@(i) uicontrol(imagePanel,'Style','text',...
        'Position',[200+30*i hPosition1 20 20],...
        'Tag',['text_channel' num2str(i)],'String',i),...
        1:numel(userData.MO.channels_));
    imagePanelSize = getPanelSize(imagePanel);
else
    imagePanel=-1;
    imagePanelSize= [0 0];
end


%% Overlay panel creation
if ~isempty(overlayProc)
    overlayPanel = uipanel(mainFig,'Position',[1/2 0 1/2 1],...
        'Title','Overlay','BackgroundColor',get(0,'defaultUicontrolBackgroundColor'),...
        'Units','pixels','Tag','uipanel_overlay');
    
    % First create overlay option (vectorField)
    hPosition2=10;
    if isempty(userData.MO.pixelSize_) || isempty(userData.MO.timeInterval_),
        scaleBarStatus = 'off';
    else
        scaleBarStatus = 'on';
    end
    uicontrol(overlayPanel,'Style','edit','Position',[30 hPosition2 50 20],...
        'String','1000','BackgroundColor','white','Tag','edit_vectorFieldScaleBar',...
        'Enable',scaleBarStatus,...
        'Callback',@(h,event) setScaleBar(guidata(h),'vectorFieldScaleBar'));
    uicontrol(overlayPanel,'Style','text','Position',[85 hPosition2-2 70 20],...
        'String','nm/min','HorizontalAlignment','left');
    uicontrol(overlayPanel,'Style','checkbox',...
        'Position',[150 hPosition2 100 20],'Tag','checkbox_vectorFieldScaleBarLabel',...
        'String',' Show label','HorizontalAlignment','left',...
        'Enable',scaleBarStatus,...
        'Callback',@(h,event) setScaleBar(guidata(h),'vectorFieldScaleBar'));
    
    hPosition2=hPosition2+30;
    uicontrol(overlayPanel,'Style','checkbox',...
        'Position',[20 hPosition2 100 20],'Tag','checkbox_vectorFieldScaleBar',...
        'String',' Scalebar','HorizontalAlignment','left',...
        'Enable',scaleBarStatus,...
        'Callback',@(h,event) setScaleBar(guidata(h),'vectorFieldScaleBar'));
    uicontrol(overlayPanel,'Style','popupmenu','Position',[130 hPosition2 120 20],...
        'String',{'NorthEast', 'SouthEast', 'SouthWest', 'NorthWest'},'Value',3,...
        'Tag','popupmenu_vectorFieldScaleBarLocation','Enable',scaleBarStatus,...
        'Callback',@(h,event) setScaleBar(guidata(h),'vectorFieldScaleBar'));
    
    hPosition2=hPosition2+30;
    uicontrol(overlayPanel,'Style','text',...
        'Position',[20 hPosition2 100 20],'Tag','text_vectorFieldScale',...
        'String',' Display scale','HorizontalAlignment','left');
    uicontrol(overlayPanel,'Style','edit','Position',[120 hPosition2 50 20],...
        'String','1','BackgroundColor','white','Tag','edit_vectorFieldScale',...
        'Callback',@(h,event) redrawOverlays(guidata(h)));
    
    hPosition2=hPosition2+20;
    uicontrol(overlayPanel,'Style','text',...
        'Position',[10 hPosition2 200 20],'Tag','text_vectorFieldOptions',...
        'String','Vector field options','HorizontalAlignment','left','FontWeight','bold');
    
    % Create controls for selecting overlays
    hPosition2=hPosition2+50;
    nProc = numel(overlayProc);
    for iProc=nProc:-1:1;
        output=overlayProc{iProc}.getDrawableOutput;
        
        % Create checkboxes for movie overlays
        validOutput = find(strcmp({output.type},'movieOverlay'));
        for iOutput=validOutput(end:-1:1)
            createMovieBox(overlayPanel,overlayProcId(iProc),iOutput,hPosition2,output(iOutput).name,...
                'Callback',@(h,event) redrawOverlay(h,guidata(h)));
            hPosition2=hPosition2+20;
        end
        
        % Create checkboxes for channel-specific overlays
        validOutput = find(strcmp({output.type},'overlay'));
        for iOutput=validOutput(end:-1:1)
            validChan = overlayProc{iProc}.checkChannelOutput;
            createOutputText(overlayPanel,overlayProcId(iProc),iOutput,hPosition2,output(iOutput).name);
            arrayfun(@(x) createChannelBox(overlayPanel,overlayProcId(iProc),iOutput,x,hPosition2,...
                'Callback',@(h,event) redrawOverlay(h,guidata(h))),find(validChan));
            hPosition2=hPosition2+20;
        end
        createProcText(overlayPanel,overlayProcId(iProc),iOutput,hPosition2,overlayProc{iProc}.getName);
        hPosition2=hPosition2+20;
    end
    
    if ~isempty(overlayProc)
        uicontrol(overlayPanel,'Style','text','Position',[120 hPosition2 100 20],...
            'Tag','text_channels','String','Channels');
        arrayfun(@(i) uicontrol(overlayPanel,'Style','text',...
            'Position',[200+30*i hPosition2 20 20],...
            'Tag',['text_channel' num2str(i)],'String',i),...
            1:numel(userData.MO.channels_));
    end
    overlayPanelSize = getPanelSize(overlayPanel);
else
    overlayPanel=-1;
    overlayPanelSize= [0 0];
end
%% Add additional panel for independent graphs
if ~isempty(graphProc) 
    graphPanel = uipanel(mainFig,'Position',[0 0 1 1],...
        'Title','Graph','BackgroundColor',get(0,'defaultUicontrolBackgroundColor'),...
        'Units','pixels','Tag','uipanel_graph');
    hPosition3=10;
    
    % Create controls for selecting all other graphs
    nProc = numel(graphProc);
    for iProc=nProc:-1:1;
        output=graphProc{iProc}.getDrawableOutput;
        if isa(graphProc{iProc},'CorrelationCalculationProcess');

            validOutput = find(strcmp({output.type},'correlationGraph'));
            % Create set of boxes for correlation graphs (input/input)
            for iOutput=validOutput(end:-1:1)
                input=graphProc{iProc}.getInput;
                nInput=numel(input);
                for iInput=nInput:-1:1
                    createOutputText(graphPanel,graphProcId(iProc),iInput,hPosition3,input(iInput).name);
                    for jInput=1:iInput
                        createInputInputBox(graphPanel,graphProcId(iProc),iOutput,iInput,jInput,hPosition3,...
                            'Callback',@(h,event) redrawCorrelationGraph(h,guidata(h)));
                    end
                    hPosition3=hPosition3+20;
                end
                createProcText(graphPanel,graphProcId(iProc),iInput,hPosition3,output(iOutput).name);
                hPosition3=hPosition3+20;
            end
        elseif isa(graphProc{iProc},'EventAlignmentProcess');
            % Create set of boxes for correlation graphs (input/input)
            input=graphProc{iProc}.getInput;
            nInput=numel(input);
            validOutput = find(strcmp({output.type},'graph'));
            for iOutput=validOutput(end:-1:1)
                for iInput=nInput:-1:1                       
                    createInputBox(graphPanel,graphProcId(iProc),iOutput,iInput,hPosition3,...
                        input(iInput).name,'Callback',@(h,event) redrawEventGraph(h,guidata(h)));                  
                    hPosition3=hPosition3+20;
                end
                createProcText(graphPanel,graphProcId(iProc),iInput,hPosition3,output(iOutput).name);
                hPosition3=hPosition3+20;
            end
            
        else   
            % Create boxes for movie -specific graphs
            validOutput = find(strcmp({output.type},'movieGraph'));
            for iOutput=validOutput(end:-1:1)
                createMovieBox(graphPanel,graphProcId(iProc),iOutput,hPosition3,...
                    output(iOutput).name,'Callback',@(h,event) redrawGraph(h,guidata(h)));
                hPosition3=hPosition3+20;
            end
            
            % Create boxes for channel-specific graphs
            validOutput = find(strcmp({output.type},'graph'));
            for iOutput=validOutput(end:-1:1)
                validChan = graphProc{iProc}.checkChannelOutput();
                createOutputText(graphPanel,graphProcId(iProc),iOutput,hPosition3,output(iOutput).name);
                arrayfun(@(x) createChannelBox(graphPanel,graphProcId(iProc),iOutput,x,hPosition3,...
                    'Callback',@(h,event) redrawGraph(h,guidata(h))),find(validChan));
                hPosition3=hPosition3+20;
            end

            createProcText(graphPanel,graphProcId(iProc),iOutput,hPosition3,graphProc{iProc}.getName);
            hPosition3=hPosition3+20;
        end
       
    end
    
    if ~isempty(graphProc) && isa(userData.MO,'MovieData')
        uicontrol(graphPanel,'Style','text','Position',[120 hPosition3 100 20],...
            'Tag','text_channels','String','Channels');
        arrayfun(@(i) uicontrol(graphPanel,'Style','text',...
            'Position',[200+30*i hPosition3 20 20],...
            'Tag',['text_channel' num2str(i)],'String',i),...
            1:numel(userData.MO.channels_));
    end
    graphPanelSize = getPanelSize(graphPanel);
else
    graphPanel=-1;
    graphPanelSize= [0 0];
end


%% Get image/overlay panel size and resize them
panelsLength = max(500,imagePanelSize(1)+overlayPanelSize(1)+graphPanelSize(1));
panelsHeight = max([imagePanelSize(2),overlayPanelSize(2),graphPanelSize(2)]);

% Resize panel
if ishandle(imagePanel)
    set(imagePanel,'Position',[10 panelsHeight-imagePanelSize(2)+10 ...
        imagePanelSize(1) imagePanelSize(2)],...
        'SelectionChangeFcn',@(h,event) redrawImage(guidata(h)))
end
if ishandle(overlayPanel)
    set(overlayPanel,'Position',[imagePanelSize(1)+10 panelsHeight-overlayPanelSize(2)+10 ...
        overlayPanelSize(1) overlayPanelSize(2)]);
end
if ishandle(graphPanel)
    set(graphPanel,'Position',[imagePanelSize(1)+overlayPanelSize(1)+10 ...
        panelsHeight-graphPanelSize(2)+10 ...
        graphPanelSize(1) graphPanelSize(2)])
end

%% Create movie panel
moviePanel = uipanel(mainFig,...
    'Title','','BackgroundColor',get(0,'defaultUicontrolBackgroundColor'),...
    'Units','pixels','Tag','uipanel_movie','BorderType','none');

% Create control button for exporting figures and movie (cf Francois' GUI)
hPosition=10;

handles.movieButton = uicontrol(moviePanel, 'Style', 'pushbutton', ...
    'String', 'Make movie',...
    'Position', [150 hPosition 100 20],...
    'Callback', @(h,event) makeMovie(h,guidata(h)));

if isa(userData.MO,'MovieData')
    % Create controls for scrollling through the movie
    hPosition = hPosition+30;
    uicontrol(moviePanel,'Style','text','Position',[10 hPosition 50 15],...
        'String','Frame','Tag','text_frame','HorizontalAlignment','left');
    uicontrol(moviePanel,'Style','edit','Position',[70 hPosition 30 20],...
        'String','1','Tag','edit_frame','BackgroundColor','white',...
        'HorizontalAlignment','left',...
        'Callback',@(h,event) redrawScene(h,guidata(h)));
    uicontrol(moviePanel,'Style','text','Position',[100 hPosition 40 15],...
        'HorizontalAlignment','left',...
        'String',['/' num2str(userData.MO.nFrames_)],'Tag','text_frameMax');
    
    uicontrol(moviePanel,'Style','slider',...
        'Position',[150 hPosition panelsLength-160 20],...
        'Value',1,'Min',1,'Max',userData.MO.nFrames_,...
        'SliderStep',[1/double(userData.MO.nFrames_)  5/double(userData.MO.nFrames_)],...
        'Tag','slider_frame','BackgroundColor','white',...
        'Callback',@(h,event) redrawScene(h,guidata(h)));
end
% Create movie location edit box
hPosition = hPosition+30;
uicontrol(moviePanel,'Style','text','Position',[10 hPosition 40 20],...
    'String','Movie','Tag','text_movie');
[~,moviePath] = fileparts(ip.Results.MO.getPath);
if isa(ip.Results.MO,'MovieList')
    [~,allPaths] = cellfun(@(x) fileparts(x.getPath),userData.ML.movies_,'UniformOutput',false);
    movieIndex=0:numel(allPaths);
    uicontrol(moviePanel,'Style','popupmenu','Position',[60 hPosition panelsLength-70 20],...
        'String',vertcat(moviePath,allPaths'),'UserData',movieIndex,...
        'Value',find(userData.movieIndex==movieIndex),...
        'HorizontalAlignment','left','BackgroundColor','white','Tag','popup_movie',...
        'Callback',@(h,event) switchMovie(h,guidata(h)));
    
else
    uicontrol(moviePanel,'Style','edit','Position',[60 hPosition panelsLength-70 20],...
        'String',moviePath,...
        'HorizontalAlignment','left','BackgroundColor','white','Tag','edit_movie');

end
% Add copyrigth
hPosition = hPosition+30;
uicontrol(moviePanel,'Style','text','Position',[10 hPosition panelsLength 20],...
    'String',userfcn_softwareConfig(),'Tag','text_copyright',...
    'HorizontalAlignment','left');

% Get overlay panel size
moviePanelSize = getPanelSize(moviePanel);
moviePanelHeight =moviePanelSize(2);
set(moviePanel,'Position',[10 panelsHeight+10 panelsLength moviePanelHeight]);

%% Resize panels and figure
sz=get(0,'ScreenSize');
figWidth = panelsLength+20;
figHeight = panelsHeight+moviePanelHeight;
set(mainFig,'Position',[sz(3)/50 (sz(4)-figHeight)/2 figWidth figHeight]);


% Update handles structure and attach it to the main figure
handles = guihandles(mainFig);
guidata(handles.figure1, handles);

% Set the figure handle to -1 by default
userData.drawFig=-1;
set(handles.figure1,'UserData',userData);

%% Set up default parameters
% Auto check input process
for i=intersect(procId,validProcId)
    h=findobj(mainFig,'-regexp','Tag',['(\w)_process' num2str(i)  '_output1.*'],...
        '-not','Style','text');
    set(h,'Value',1);
    for j=find(arrayfun(@(x)isequal(get(x,'Parent'),graphPanel),h))'
        callbackFcn = get(h(j),'Callback');
        callbackFcn(h(j),[]);
    end
end

% Update the image and overlays
if isa(userData.MO,'MovieData'), redrawScene(handles.figure1, handles); end

function switchMovie(hObject,handles)
userData=get(handles.figure1,'UserData');
props=get(hObject,{'UserData','Value'});
if isequal(props{1}(props{2}), userData.movieIndex),return;end
movieViewer(userData.ML,userData.procId,'movieIndex',props{1}(props{2}));

function size = getPanelSize(hPanel)

a=get(get(hPanel,'Children'),'Position');
P=vertcat(a{:});
size = [max(P(:,1)+P(:,3))+10 max(P(:,2)+P(:,4))+20];


function makeMovie(hObject,handles)

userData = get(handles.figure1, 'UserData');
nFrames = userData.MO.nFrames_;

fmt = ['%0' num2str(ceil(log10(nFrames))) 'd'];
frameName = @(frame) ['frame' num2str(frame, fmt) '.png'];
fpath = [userData.MO.outputDirectory_ filesep 'Frames'];
mkClrDir(fpath);
fprintf('Generating movie frames:     ');
resolution = ['-r' num2str(5*72)];
for iFrame=1:nFrames
    set(handles.slider_frame, 'Value',iFrame);
    redrawScene(hObject, handles);
    drawnow;
    print(userData.drawFig, '-dpng', '-loose', resolution, fullfile(fpath,frameName(iFrame)));
    fprintf('\b\b\b\b%3d%%', round(100*iFrame/(nFrames)));
end
fprintf('\n');

% Generate movie
mpath = [userData.MO.outputDirectory_ filesep 'Movie'];
mkClrDir(mpath);
fprintf('Generating movie... ');
fr = num2str(15);
cmd = ['ffmpeg -y -r ' fr ' -i ' fpath 'frame' fmt '.png' ' -r ' fr ' -b 50000k -bt 20000k ' mpath 'movie.mp4 > /dev/null 2>&1' ];
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
frameNumber = min(max(frameNumber,1),userData.MO.nFrames_);

% Set the slider and editboxes values
set(handles.edit_frame,'String',frameNumber);
set(handles.slider_frame,'Value',frameNumber);

% Update the image and overlays
redrawImage(handles);
redrawOverlays(handles);

function h= getFigure(handles,figName)

h = findobj(0,'-regexp','Name',['^' figName '$']);
if ~isempty(h), figure(h); return; end

%Create a figure
if strcmp(figName,'Movie')
    userData = get(handles.figure1,'UserData');
    sz=get(0,'ScreenSize');
    nx=userData.MO.imSize_(2);
    ny=userData.MO.imSize_(1);
    h = figure('Position',[sz(3)*.2 sz(4)*.2 nx ny],...
        'Name',figName,'NumberTitle','off','Tag','viewerFig');
    
    % figure options for movie export
    iptsetpref('ImshowBorder','tight');
    set(h, 'InvertHardcopy', 'off');
    set(h, 'PaperUnits', 'Points');
    set(h, 'PaperSize', [nx ny]);
    set(h, 'PaperPosition', [0 0 nx ny]); % very important
    set(h, 'PaperPositionMode', 'auto');
    % set(h,'DefaultLineLineSmoothing','on');
    % set(h,'DefaultPatchLineSmoothing','on');
    
    axes('Parent',h,'XLim',[0 userData.MO.imSize_(2)],...
        'YLim',[0 userData.MO.imSize_(1)],'Position',[0.05 0.05 .9 .9]);
    userData.drawFig=h;
    set(handles.figure1,'UserData',userData);
else
    h = figure('Name',figName,'NumberTitle','off','Tag','viewerFig');
end


function redrawChannel(hObject,handles)

% Callback for channels checkboxes to avoid 0 or more than 4 channels
channelBoxes = findobj(handles.figure1,'-regexp','Tag','checkbox_channel*');
nChan=numel(find(arrayfun(@(x)get(x,'Value'),channelBoxes)));
if nChan==0, set(hObject,'Value',1); elseif nChan>3, set(hObject,'Value',0); end

redrawImage(handles)

function setScaleBar(handles,type)
% Remove existing scalebar of given type
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
    width = scale *1000/userData.MO.pixelSize_;
    label = [num2str(scale) ' \mum'];
else
    displayScale = str2double(get(handles.edit_vectorFieldScale,'String'));
    width = scale*displayScale/(userData.MO.pixelSize_/userData.MO.timeInterval_*60);
    label= [num2str(scale) ' nm/min'];
end
if ~get(handles.(['checkbox_' type 'Label']),'Value'), label=''; end
props=get(handles.(['popupmenu_' type 'Location']),{'String','Value'});
location=props{1}{props{2}};
hScaleBar = plotScaleBar(width,'Label',label,'Location',location);
set(hScaleBar,'Tag',type);

function setTimeStamp(handles)
% Remove existing timestamp of given type
h=findobj('Tag','timeStamp');
if ~isempty(h), delete(h); end

% If checked, adds a new scalebar using the width as a label input
userData=get(handles.figure1,'UserData');
if ~get(handles.checkbox_timeStamp,'Value') || ~ishandle(userData.drawFig),
    return
end
figure(userData.drawFig)
frameNr=get(handles.slider_frame,'Value');
width = userData.MO.imSize_(2)/20;
time= (frameNr-1)*userData.MO.timeInterval_;
p=sec2struct(time);
props=get(handles.popupmenu_timeStampLocation,{'String','Value'});
location=props{1}{props{2}};
hTimeStamp = plotScaleBar(width,'Label',p.str,'Location',location);
set(hTimeStamp,'Tag','timeStamp');
delete(hTimeStamp(1))

function setCLim(handles)
userData=get(handles.figure1,'UserData');
imageTag = get(get(handles.uipanel_image,'SelectedObject'),'Tag');

clim=[str2double(get(handles.edit_cmin,'String')) ...
    str2double(get(handles.edit_cmax,'String'))];
redrawImage(handles,'CLim',clim)


function setColormap(handles)
allCmap=get(handles.popupmenu_colormap,'String');
selectedCmap = get(handles.popupmenu_colormap,'Value');
redrawImage(handles,'Colormap',allCmap{selectedCmap})

function setColorbar(handles)
cbar=get(handles.checkbox_colorbar,'Value');
if cbar, cbarStatus='on'; else cbarStatus='off'; end 
redrawImage(handles,'Colorbar',cbarStatus)

function redrawImage(handles,varargin)
frameNr=get(handles.slider_frame,'Value');
imageTag = get(get(handles.uipanel_image,'SelectedObject'),'Tag');

% Get the figure handle
getFigure(handles,'Movie');
userData=get(handles.figure1,'UserData');

% Use corresponding method depending if input is channel or process output
channelBoxes = findobj(handles.figure1,'-regexp','Tag','checkbox_channel*');
[~,index]=sort(arrayfun(@(x) get(x,'Tag'),channelBoxes,'UniformOutput',false));
channelBoxes =channelBoxes(index);
if strcmp(imageTag,'radiobutton_channels')
    set(channelBoxes,'Enable','on');
    chanList=find(arrayfun(@(x)get(x,'Value'),channelBoxes));
    userData.MO.channels_(chanList).draw(frameNr,varargin{:});
    displayMethod = userData.MO.channels_(chanList(1)).displayMethod_;
else
    set(channelBoxes,'Enable','off');
    % Retrieve the id, process nr and channel nr of the selected imageProc
    tokens = regexp(imageTag,'radiobutton_process(\d+)_output(\d+)_channel(\d+)','tokens');
    procId=str2double(tokens{1}{1});
    outputList = userData.MO.processes_{procId}.getDrawableOutput;
    iOutput = str2double(tokens{1}{2});
    output = outputList(iOutput).var;
    iChan = str2double(tokens{1}{3});
    userData.MO.processes_{procId}.draw(iChan,frameNr,'output',output,varargin{:});
    displayMethod = userData.MO.processes_{procId}.displayMethod_{iOutput,iChan};
end


% Set the color limits properties
clim=displayMethod.CLim;
if isempty(clim)
    userData = get(handles.figure1,'UserData');
    hAxes=findobj(userData.drawFig,'Type','axes','-not','Tag','Colorbar');
    clim=get(hAxes,'Clim');
end
set(handles.edit_cmin,'Enable','on','String',clim(1));
set(handles.edit_cmax,'Enable','on','String',clim(2));
    
% Set the colorbar properties
cbar=displayMethod.Colorbar;
set(handles.checkbox_colorbar,'Value',strcmpi(cbar,'on'));

% Set the colormap properties
cmap=displayMethod.Colormap;
allCmap=get(handles.popupmenu_colormap,'String');
set(handles.popupmenu_colormap,'Value',find(strcmpi(cmap,allCmap)));

% Reset the scaleBar
setScaleBar(handles,'imageScaleBar');
setTimeStamp(handles);

function redrawOverlays(handles)
overlayBoxes = findobj(handles.uipanel_overlay,'-regexp','Tag','checkbox_process*');
checkedBoxes = logical(arrayfun(@(x) get(x,'Value'),overlayBoxes));
overlayTags=arrayfun(@(x) get(x,'Tag'),overlayBoxes(checkedBoxes),...
    'UniformOutput',false);
for i=1:numel(overlayTags),
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

% Get figure handle or recreate figure
if ishandle(userData.drawFig),
    figure(userData.drawFig);
else
    redrawScene(hObject, handles); return;
end
% Retrieve the id, process nr and channel nr of the selected imageProc
tokens = regexp(overlayTag,'^checkbox_process(\d+)_output(\d+)','tokens');
procId=str2double(tokens{1}{1});
outputList = userData.MO.processes_{procId}.getDrawableOutput;
iOutput = str2double(tokens{1}{2});
output = outputList(iOutput).var;

% Discriminate between channel-specific processes annd movie processes
tokens = regexp(overlayTag,'_channel(\d+)$','tokens');
if ~isempty(tokens)
    iChan = str2double(tokens{1}{1});
    inputArgs={iChan,frameNr};
    graphicTag =[userData.MO.processes_{procId}.getName '_channel'...
        num2str(iChan) '_output' num2str(iOutput)];
else
    inputArgs={frameNr};
    graphicTag = [userData.MO.processes_{procId}.getName '_output' num2str(iOutput)];
    
end

% Draw or delete the overlay depending on the checkbox value
if get(hObject,'Value')
    userData.MO.processes_{procId}.draw(inputArgs{:},'output',output,...
        'vectorScale',str2double(get(handles.edit_vectorFieldScale,'String')));
else
    h=findobj('Tag',graphicTag);
    if ~isempty(h), delete(h); end
end

function redrawGraph(hObject,handles)
overlayTag = get(hObject,'Tag');
userData=get(handles.figure1,'UserData');

% Retrieve the id, process nr and channel nr of the selected graphProc
tokens = regexp(overlayTag,'^checkbox_process(\d+)_output(\d+)','tokens');
procId=str2double(tokens{1}{1});
outputList = userData.MO.processes_{procId}.getDrawableOutput;
iOutput = str2double(tokens{1}{2});
output = outputList(iOutput).var;

% Discriminate between channel-specific and movie processes
tokens = regexp(overlayTag,'_channel(\d+)$','tokens');
if ~isempty(tokens)
    iChan = str2double(tokens{1}{1});
    inputArgs={iChan};
    figName = [outputList(iOutput).name ' - Channel ' num2str(iChan)];
else
    inputArgs={};
    figName = outputList(iOutput).name;
end

% Draw or delete the graph figure depending on the checkbox value
if get(hObject,'Value')
    h = getFigure(handles,figName);
    userData.MO.processes_{procId}.draw(inputArgs{:},'output',output,...
        'vectorScale',str2double(get(handles.edit_vectorFieldScale,'String')));
    set(h,'DeleteFcn',@(h,event)closeGraphFigure(hObject));
else
    h=findobj(0,'-regexp','Name',['^' figName '$']);
    if ~isempty(h), delete(h); end
end



function redrawCorrelationGraph(hObject,handles)
overlayTag = get(hObject,'Tag');
userData=get(handles.figure1,'UserData');

% Retrieve the id, process nr and channel nr of the selected graphProc
tokens = regexp(overlayTag,'^checkbox_process(\d+)_output(\d+)_input(\d+)_input(\d+)','tokens');
procId=str2double(tokens{1}{1});
outputList = userData.MO.processes_{procId}.getDrawableOutput;
input = userData.MO.processes_{procId}.getInput;
iOutput = str2double(tokens{1}{2});
iInput1 = str2double(tokens{1}{3});
iInput2 = str2double(tokens{1}{4});
output = outputList(iOutput).var;

if iInput1==iInput2
    figName = [input(iInput1).name ' autocorrelation'];
else
    figName = [input(iInput1).name ' - ' input(iInput2).name ' cross-correlation'];
end


% Draw or delete the graph figure depending on the checkbox value
if get(hObject,'Value')
    h = getFigure(handles,figName);
    userData.MO.processes_{procId}.draw(iInput1,iInput2,'output',output);
    set(h,'DeleteFcn',@(h,event)closeGraphFigure(hObject));
else
    h=findobj(0,'-regexp','Name',['^' figName '$']);
    if ~isempty(h), delete(h); end
end

function redrawEventGraph(hObject,handles)
overlayTag = get(hObject,'Tag');
userData=get(handles.figure1,'UserData');

% Retrieve the id, process nr and channel nr of the selected graphProc
tokens = regexp(overlayTag,'^checkbox_process(\d+)_output(\d+)_input(\d+)','tokens');
procId=str2double(tokens{1}{1});
outputList = userData.MO.processes_{procId}.getDrawableOutput;
input = userData.MO.processes_{procId}.getInput;
iOutput = str2double(tokens{1}{2});
iInput1 = str2double(tokens{1}{3});
output = outputList(iOutput).var;


figName = ['Aligned ' input(iInput1).name];

% Draw or delete the graph figure depending on the checkbox value
if get(hObject,'Value')
    h = getFigure(handles,figName);
    userData.MO.processes_{procId}.draw(iInput1);
    set(h,'DeleteFcn',@(h,event)closeGraphFigure(hObject));
else
    h=findobj(0,'-regexp','Name',['^' figName '$']);
    if ~isempty(h), delete(h); end
end

function closeGraphFigure(hObject)
set(hObject,'Value',0);

function deleteViewer()

h = findobj(0,'-regexp','Tag','viewerFig');
if ~isempty(h), delete(h); end
