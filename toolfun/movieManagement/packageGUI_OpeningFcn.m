function packageGUI_OpeningFcn(hObject,eventdata,handles,packageName,varargin)
% Callback called at the opening of packageGUI
%
% packageGUI(MD)   MD: MovieData object
%
% Useful tools
%
% User Data:
%
%       userData.MD - array of MovieData object
%       userData.package - array of package (same length with userData.MD)
%       userData.crtPackage - the package of current MD
%       userData.id - the id of current MD on board
%
%       userData.dependM - dependency matrix
%       userdata.statusM - GUI status matrix
%       userData.optProcID - optional process ID
%       userData.applytoall - array of boolean for batch movie set up
%
%       userData.passIconData - pass icon image data
%       userData.errorIconData - error icon image data
%       userData.warnIconData - warning icon image data
%       userData.questIconData - help icon image data
%       userData.colormap - color map
%
%       userData.setFig - array of handles of (multiple) setting figures (may not exist)
%       userData.resultFig - array of handles of (multiple) result figures (may not exist)
%       userData.packageHelpFig - handle of (single) help figure (may not exist)
%       userData.iconHelpFig - handle of (single) help figures (may not exist)
%       userData.processHelpFig - handle of (multiple) help figures (may not exist) 
%       
%
% NOTE:
%   
%   userData.statusM - 1 x m stucture array, m is the number of Movie Data 
%                      this user data is used to save the status of movies
%                      when GUI is switching between different movie(s)
%                   
%   	fields: IconType - the type of status icons, 'pass', 'warn', 'error'
%               Msg - the message displayed when clicking status icons
%               Checked - 1 x n logical array, n is the number of processes
%                         used to save value of check box of each process
%               Visited - logical true or false, if the movie has been
%                         loaded to GUI before 
% Sebastien Besson 5/2011

ip = inputParser;
ip.addRequired('hObject',@ishandle);
ip.addRequired('eventdata',@(x) isstruct(x) || isempty(x));
ip.addRequired('handles',@isstruct);
ip.addRequired('packageName',@(x) isa(x,'char') || isa(x,'function_handle'));
ip.addOptional('MD',[],@(x) isa(x,'MovieData'));
ip.parse(hObject,eventdata,handles,packageName,varargin{:});
MD=ip.Results.MD;

if isa(packageName,'function_handle')
    packageHandle=packageName;
    packageName=func2str(packageName);
elseif isa(packageName,'char')
    packageHandle=str2func(packageName);
end

assert(any(strcmp(superclasses(packageName),'Package')),sprintf('%s is not a valid Package',packageName));
      
handles.output = hObject;
userData = get(handles.figure1,'UserData');
userData.packageName = packageName;


% Call package GUI error

[copyright openHelpFile] = userfcn_softwareConfig(handles);
set(handles.text_copyright, 'String', copyright);

%If package GUI supplied without argument, saves a boolean which will be
%read by packageNameGUI_OutputFcn
if isempty(MD)
    userData.startMovieSelectorGUI=true;
    set(handles.figure1,'UserData',userData);
    guidata(hObject, handles);
    return
end

% ----------------------------- Load MovieData ----------------------------

% MD = varargin{2};
nMovies = numel(MD);
packageIndx = cell(1, nMovies);

% I. Before loading MovieData, firstly check if the current package exists
for i = 1:nMovies
    % Check for existing packages and create them if false
    packageIndx{i} = find(cellfun(@(x) isa(x,packageName),MD(i).packages_),1);
    if packageIndx{i}
        userData.package(i) = MD(i).packages_{packageIndx{i}};
    else
        MD(i).addPackage(packageHandle(MD(i), MD(i).outputDirectory_));
        userData.package(i) = MD(i).packages_{end};
        % Sanity check to check basic dependencies are satisfied
        try
            userData.package(i).sanityCheck(true,'all');
        catch ME
            errordlg(ME.message);
            userData.startMovieSelectorGUI=true;
            set(handles.figure1,'UserData',userData);
            guidata(hObject, handles);
            return
        end
    end
end

% ------------- Check if existing processes can be recycled ---------------

existProcess = cell(1, nMovies);
processClassNames = userData.package(1).processClassNames_;

% Multiple movies loop
for i = 1:nMovies

    if isempty(packageIndx{i}) && ~isempty(MD(i).processes_)
    
        classname = cellfun(@(z)class(z), MD(i).processes_, 'UniformOutput', false);

        existProcessForm = cellfun(@(z)strcmp(z, classname), processClassNames, 'uniformoutput', false );
        existProcessId = find(cellfun(@(z)any(z), existProcessForm));
    
        if ~isempty (existProcessId)
            % Get recycle processes 
        
            for j = existProcessId
                existProcess{i} = horzcat(existProcess{i},  MD(i).processes_(existProcessForm{j}) );
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

% Initialize userdata
userData.id = 1;
userData.crtPackage = userData.package(userData.id);
userData.MD = MD;
userData.dependM = userData.package(userData.id).getDependencyMatrix;
userData.optProcID =userData.package(userData.id).getOptionalProcessId;
nProc = size(userData.dependM, 1);
userData.statusM = repmat( struct('IconType', {cell(1,nProc)}, 'Msg', {cell(1,nProc)}, 'Checked', zeros(1,nProc), 'Visited', false), 1, nMovies);

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
set(handles.figure1,'CurrentAxes',handles.axes_help);
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

% Set option depen
if length(userData.MD) == 1
    set(handles.checkbox_runall, 'Visible', 'off')
    set(handles.pushbutton_left, 'Enable', 'off')
    set(handles.pushbutton_right, 'Enable', 'off')   
    set(handles.checkbox_all, 'Visible', 'off', 'Value', 0)
    userData.applytoall=zeros(nProc,1);
else
    set(handles.checkbox_runall, 'Visible', 'on')
    userData.applytoall=ones(nProc,1);
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


end