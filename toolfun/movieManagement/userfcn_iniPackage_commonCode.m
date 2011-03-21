% Not a function, Code segment

stack = dbstack;
guiname=stack(4).name;
packageName=regexprep(guiname(1:end-3),'(\<[a-z])','${upper($1)}');

% if  strcmp(stack(2).name,'MovieData.relocateMovieData')
handles.output = hObject;
userData = get(handles.figure1,'UserData');

eval(['userData.optProcID =' packageName '.getOptionalProcessId;']);

% Call package GUI error

[copyright openHelpFile] = userfcn_softwareConfig(handles);
set(handles.text_copyright, 'String', copyright)

if nargin < 4
    handles.startMovieSelectorGUI=true;
    handles.packageName = packageName;
    guidata(hObject, handles);
    return
%     error('User-defined: Please call package control panel with a MovieData object. E.g. packageGUI(movieDataObject)');
end

% ----------------------------- Load MovieData ----------------------------

MD = varargin{1};
len = length(MD);
packageExist = zeros(1, len);
% I. Before loading MovieData, firstly check if the current package exists


for x = 1:len
    
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
        eval(['MD(x).addPackage( ' packageName '(MD(x), MD(x).outputDirectory_) )'])
        userData.package(x) = MD(x).packages_{end};
    end

end

% ------------- Check if existing processes can be recycled ---------------

existProcess = cell(1, len);
processClassNames = userData.package(1).processClassNames_;

% Multiple movies loop
for x = 1:len

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
eval(['userData.dependM =  ' packageName '.getDependencyMatrix;'])

l = size(userData.dependM, 1);
userData.statusM = repmat( struct('IconType', {cell(1,l)}, 'Msg', {cell(1,l)}, 'Checked', zeros(1,l), 'Visited', false), 1, len);



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
% set(Img,'ButtonDownFcn',@userfcn_openHelp('UTrackPackage'));
set(Img,'ButtonDownFcn',@icon_ButtonDownFcn);

if openHelpFile
    set(Img, 'UserData', struct('class', packageName))
end

% Set up process help
for i = 1:l

    eval (['axes(handles.axes_help_' num2str(i) ')'])
    Img = image(questIconData);
    set(gca, 'XLim',get(Img,'XData'),'YLim',get(Img,'YData'),...
        'visible','off','YDir','reverse');  
    set(Img,'ButtonDownFcn',@icon_ButtonDownFcn);
        
    if openHelpFile
        set(Img, 'UserData', struct('class', userData.crtPackage.processClassNames_{i}))
    end

end

% --------------------------Other GUI settings-----------------------------

% set text body
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
    set(handles.checkbox_all, 'Visible', 'off', 'Value', 0)
else
    set(handles.checkbox_runall, 'Visible', 'on')
end

% Update handles structure
set(handles.figure1,'UserData',userData);
guidata(hObject, handles);

userfcn_updateGUI(handles, 'initialize')
 
