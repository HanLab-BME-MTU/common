
function userfcn_updateGUI(handles, type)
% GUI tool function: this function is called by movie explorer when 
% switching between differenct movies. 
% 
% Input: 
%       handles - the "handles" of package GUI control panel
%       type - 'initialize': used when movie is loaded to GUI for the first time
%              'refresh': used when movie had already been loaded to GUI
%
%
% Chuangang Ren
% 08/2010

if ~any(strcmp(type, {'initialize' 'refresh'}))
    error('User-defined: Error calling userfcn_updateGUI, type should be ''initialize'' or ''refresh''')
end

userData = get(handles.figure1, 'UserData');
l = size(userData.dependM, 1);
k = zeros(1,l);

% reset GUI
userfcn_resetGUI(handles)

% Set movie data path
set(handles.edit_path, 'String', ...
    [userData.MD(userData.id).movieDataPath_ filesep userData.MD(userData.id).movieDataFileName_ ])


% ----------------------------- Initialize --------------------------------

if strcmp(type, 'initialize')

    % Package Sanity Check

    procEx = userData.crtPackage.sanityCheck(true, 'all');

    for i = 1: l

        % Return user data !!!
        set(handles.figure1, 'UserData', userData)

       if ~isempty(procEx{i})

           if strcmp(procEx{i}(1).identifier, 'lccb:set:fatal')
               userfcn_drawIcon(handles,'error',i,procEx{i}(1).message, true); % user data is retrieved, updated and submitted
           else
               userfcn_drawIcon(handles,'warn',i,procEx{i}(1).message, true); % user data is retrieved, updated and submitted
           end

       else
           if ~isempty(userData.crtPackage.processes_{i}) && ...
              userData.crtPackage.processes_{i}.success_ && ...
               ~userData.crtPackage.processes_{i}.procChanged_ && ...
               userData.crtPackage.processes_{i}.updated_

               userfcn_drawIcon(handles,'pass',i,'Current step was processed successfully', true) ; % user data is retrieved, updated and submitted

           end
       end

       % Refresh user data !!!
        userData = get(handles.figure1, 'UserData');
    end
    
% ----------------------------- Refresh -----------------------------------    
elseif strcmp(type, 'refresh')
    for i = 1: l
    
       % Draw icons
       if ~isempty(userData.statusM(userData.id).IconType{i})
            userfcn_drawIcon(handles, userData.statusM(userData.id).IconType{i}, i, userData.statusM(userData.id).Msg{i}, false);
       end

    end
end  
% -------------------------------------------------------------------------

for i = 1: l
    
    % If process is checked, check and enable the process and enable decendent
    % processes
    if userData.statusM(userData.id).Checked(i)
        k(i) = 1;
        eval([ 'set(handles.checkbox_',num2str(i),', ''Value'',1, ''Enable'', ''on'')' ])
        userfcn_lampSwitch(i, 1, handles)
    end      
    
    % Bold the Name of Existing Process 

   if ~isempty(userData.crtPackage.processes_{i})
       eval([ 'set(handles.checkbox_',num2str(i),', ''FontWeight'',''bold'')' ]);
   end

    
    % Set Up Uicontrols Enable/Disable

   % If process's sucess = 1, release the process from GUI enable/disable
   % control
   if ~isempty(userData.crtPackage.processes_{i}) && ...
      userData.crtPackage.processes_{i}.success_ 
       k(i) = 1;
       eval([ 'set(handles.pushbutton_show_',num2str(i),', ''enable'', ''on'');']);
   end
   
end

tempDependM = userData.dependM;
tempDependM(:,logical(k)) = zeros(l, nnz(k));

% Checkbox enable/disable set up
userfcn_enable(find (any(tempDependM,2)), 'off',handles);


if strcmp(type, 'initialize')
    
    userData.statusM(userData.id).Visited = true;
    set(handles.figure1, 'UserData', userData)

end
