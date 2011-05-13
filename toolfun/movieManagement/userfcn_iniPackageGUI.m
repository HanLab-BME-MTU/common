%Package GUI initialization - should be called at the opening of the
%corresponding GUI

%Reed package name 
stack = dbstack;
if ~strcmp(stack(2).name(end-9:end),'OpeningFcn'), 
    error('lccb:Package:Initialization',[mfilename ' should be called at the GUI opening'])
end
guiname=stack(4).name;

if isempty(varargin)
    error('lccb:packageGUI',sprintf('%s must be called with at least one valid argument',mfilename));
end

if isa(varargin{1},'function_handle')
    packageHandle=varargin{1};
    packageName=func2str(varargin{1});
elseif isa(varargin{1},'char')
    packageName=varargin{1};
    packageHandle=str2func(varargin{1});
else
    packageName=class(varargin{1});
    packageHandle=str2func(packageName);
end

assert(logical(exist(packageName,'class')),sprintf('%s is not a valid class',packageName));
assert(any(strcmp(superclasses(packageName),'Package')),sprintf('%s is not a valid Package class',packageName));
      
handles.output = hObject;
userData = get(handles.figure1,'UserData');
userData.packageName = packageName;

userData.optProcID =eval([userData.packageName,'.getOptionalProcessId']);

% Call package GUI error

[copyright openHelpFile] = userfcn_softwareConfig(handles);
set(handles.text_copyright, 'String', copyright);

%If package GUI supplied without argument, saves a boolean which will be
%read by packageNameGUI_OutputFcn
if nargin < 5
    userData.startMovieSelectorGUI=true;
    set(handles.figure1,'UserData',userData);
    guidata(hObject, handles);
    return
end

% ----------------------------- Load MovieData ----------------------------

MD = varargin{2};
nMovies = numel(MD);
packageExist = zeros(1, nMovies);
% I. Before loading MovieData, firstly check if the current package exists


for x = 1:nMovies
    
    packageExist(x) = false;

    for i = 1: length(MD(x).packages_)
        if isa(MD(x).packages_{i}, packageName)

            userData.package(x) = MD(x).packages_{i};
            packageExist(x) = true;
            break;
        end
    end
    
    if ~packageExist(x)
        % No same package is found.. create a new package object
        MD(x).addPackage(packageHandle(MD(x), MD(x).outputDirectory_)) 
        userData.package(x) = MD(x).packages_{end};
    end

end

% ------------- Check if existing processes can be recycled ---------------

existProcess = cell(1, nMovies);
processClassNames = userData.package(1).processClassNames_;

% Multiple movies loop
for x = 1:nMovies

    if ~packageExist(x) && ~isempty(MD(x).processes_)
    
        classname = cellfun(@(z)class(z), MD(x).processes_, 'UniformOutput', false);

        existProcessForm = cellfun(@(z)strcmp(z, classname), processClassNames, 'uniformoutput', false );
        existProcessId = find(cellfun(@(z)any(z), existProcessForm));
    
        if ~isempty (existProcessId)
            % Get recycle processes 
        
            for i = existProcessId
                existProcess{x} = horzcat(existProcess{x},  MD(x).processes_(existProcessForm{i}) );
            end

        end
    
    end
end

existProcessMovieId = find( cellfun(@(z)(~isempty(z)), existProcess));

if ~isempty(existProcessMovieId)
    % Get messages
    procMsg = [];
    for i = 1:length(existProcess{existProcessMovieId(1)})
        procMsg = [procMsg sprintf('%s Step\n', existProcess{existProcessMovieId(1)}{i}.name_)];
    end

    msg = sprintf('Record indicates that the following steps are availabe for %s package: \n\n%s\nDo you want to load and re-use these steps in %s package?', ...
                   userData.package(1).name_, procMsg, userData.package(1).name_);
                  
    % Ask user if to recycle
    user_response = questdlg(msg, 'Recycle Existing Steps',  'No', 'Yes','Yes');
    
    if strcmpi( user_response , 'Yes')

        for x = existProcessMovieId
            
            recycleProcessGUI(existProcess{x}, userData.package(x), 'mainFig', handles.figure1)
        end
    end
        
        
end


userData.id = 1;
userData.crtPackage = userData.package(userData.id);
userData.MD = MD;

% Dependency matrix is defined in Package class
% Make a copy of dependency matrix here to control enable/disable 
% machanism in package GUI
userData.dependM = eval([packageName,'.getDependencyMatrix']);

nProc = size(userData.dependM, 1);
userData.statusM = repmat( struct('IconType', {cell(1,nProc)}, 'Msg', {cell(1,nProc)}, 'Checked', zeros(1,nProc), 'Visited', false), 1, nMovies);

% Initialize the apply to all checkboxes
userData.applytoall=zeros(nProc,1);

% -----------------------Load and set up icons----------------------------

% Load icon images from dialogicons.mat
load lccbGuiIcons.mat

% Save Icon data to GUI data
userData.passIconData = passIconData;
userData.errorIconData = errorIconData;
userData.warnIconData = warnIconData;
userData.questIconData = questIconData;

% Set figure colormap
supermap(1,:) = get(hObject,'color');
set(hObject,'colormap',supermap);

userData.colormap = supermap;

% Set up package help. 
axes(handles.axes_help);
Img = image(questIconData); 
set(gca, 'XLim',get(Img,'XData'),'YLim',get(Img,'YData'),...
    'visible','off','YDir','reverse');
set(Img,'ButtonDownFcn',@icon_ButtonDownFcn);

if openHelpFile
    set(Img, 'UserData', struct('class', packageName))
end
% --------------------------Set up processes------------------------------

% List of template process uicontrols to expand
templateTag{1} = 'checkbox';
templateTag{2} = 'axes_icon';
templateTag{3} = 'pushbutton_show';
templateTag{4} = 'pushbutton_set';
templateTag{5} = 'axes_prochelp';
% templateTag{6} = 'pushbutton_clear'; To be implemented someday?
procTag=templateTag;
set(handles.figure1,'Position',...
    get(handles.figure1,'Position')+(nProc-1)*[0 0 0 40])
set(handles.panel_movie,'Position',...
    get(handles.panel_movie,'Position')+(nProc-1)*[0 40 0 0])
set(handles.panel_proc,'Position',...
    get(handles.panel_proc,'Position')+(nProc-1)*[0 0 0 40])
set(handles.text_status, 'Position',...
    get(handles.text_status,'Position')+(nProc-1)*[0 40 0 0])      

for i = 1:nProc
    for j=1:length(templateTag)
        procTag{j}=[templateTag{j} '_' num2str(i)];
        handles.(procTag{j}) = copyobj(handles.(templateTag{j}),handles.panel_proc);
        set(handles.(procTag{j}),'Tag',procTag{j},'Position',...
            get(handles.(templateTag{j}),'Position')+(nProc-i)*[0 40 0 0]);
    end
  
    processName=userData.crtPackage.processClassNames_{i};
    checkboxString = [' Step ' num2str(i) ':' regexprep(processName,'([A-Z])',' $1')];
    set(handles.(procTag{1}),'String',checkboxString)
    
    set(handles.figure1,'CurrentAxes',handles.(procTag{5}));
    Img = image(questIconData);
    set(gca, 'XLim',get(Img,'XData'),'YLim',get(Img,'YData'),...
        'visible','off','YDir','reverse');  
    set(Img,'ButtonDownFcn',@icon_ButtonDownFcn);
        
    if openHelpFile
        set(Img, 'UserData', struct('class', processName))
    end
end

cellfun(@(x)delete(handles.(x)),templateTag)
handles = rmfield(handles,templateTag);

optTag = 'text_optional';
for i = userData.optProcID
    procOptTag=[optTag '_' num2str(i)];
    handles.(procOptTag) = copyobj(handles.(optTag),handles.panel_proc);
    set(handles.(procOptTag),'Tag',procOptTag,'Position',...
        get(handles.(optTag),'Position')+(nProc-i)*[0 40 0 0]);
end

delete(handles.(optTag));
handles = rmfield(handles,optTag);


% --------------------------Create tools menu-----------------------------

if ~isempty(userData.crtPackage.tools_)
    handles.menu_tools = uimenu(handles.figure1,'Label','Tools','Position',2);
    for i=1:length(userData.crtPackage.tools_)
        toolMenuTag=['menu_tools_' num2str(i)];
        handles.(toolMenuTag) = uimenu(handles.menu_tools,...
            'Label',userData.crtPackage.tools_(i).name,...
            'Callback',@menu_tools_Callback,'Tag',toolMenuTag);
    end
end

% --------------------------Other GUI settings-----------------------------

% set titles
set(handles.figure1, 'Name',['Control Panel - ' userData.crtPackage.name_]);
set(handles.text_body1, 'string',[userData.crtPackage.name_ ' Package']);

% Set movie explorer
msg = {};
for i = 1: length(userData.MD)
    msg = cat(2, msg, {sprintf('  Movie %d of %d', i, length(userData.MD))});
end
set(handles.popupmenu_movie, 'String', msg, 'Value', userData.id);

% If only one movie loaded

if length(userData.MD) == 1
    set(handles.checkbox_runall, 'Visible', 'off')
    set(handles.pushbutton_left, 'Enable', 'off')
    set(handles.pushbutton_right, 'Enable', 'off')   
%     set(handles.checkbox_all, 'Visible', 'off', 'Value', 0)
else
    set(handles.checkbox_runall, 'Visible', 'on')
end


% Set web links in menu
set(handles.menu_about_gpl,'UserData','http://www.gnu.org/licenses/gpl.html')
set(handles.menu_about_lccb,'UserData','http://lccb.hms.harvard.edu/')
set(handles.menu_about_lccbsoftware,'UserData','http://lccb.hms.harvard.edu/software.html')
 
% Update handles structure
set(handles.figure1,'UserData',userData);
guidata(hObject, handles);
set(Img,'ButtonDownFcn',@icon_ButtonDownFcn);

userfcn_updateGUI(handles, 'initialize')


