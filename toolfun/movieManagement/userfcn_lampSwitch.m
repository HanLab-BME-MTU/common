function userfcn_lampSwitch(index, value, handles)
% GUI tool function: control the enable/disable value of uicontrols when
% user check/unchecked the checkboxes of processes. The enable/disable 
% value of uicontrols depends on package's dependency matrix - dependM
%
% Input: 
%   index - the index of current checkbox
%   value - 1: checked   0: unchecked
%   handles - the "handles" of package control panel movie 
%
%
% Chuangang Ren
% 08/2010

userData = get(handles.figure1, 'UserData');
M = userData.dependM;


if ~any(M(:,index))
   % if no follower exists, return.
        return;
else
    childProcesses = find(M(:,index));
    switch value
        % Checkbox is selected
        case 1
            for i = 1: length(childProcesses)
               requiredParentProcesses = find(M(childProcesses(i),:)==1);
               for j = 1: length(requiredParentProcesses)
                   if get(handles.(['checkbox_' num2str(requiredParentProcesses(j))]),'Value') ||...
                           ( ~isempty(userData.crtPackage.processes_{requiredParentProcesses(j)}) && ...
                           userData.crtPackage.processes_{requiredParentProcesses(j)}.success_ ) ||...
                           j == index
                       
                       k = true; % ok
                       
                   else
                       k = false; % not ok
                       break
                   end
               end
               if ~k
                   continue;
               end
               % The following code will probably not be executed
               % Leave it here just in case design is changed
               % ------------------------------------------ %
               if get(handles.(['checkbox_' num2str(childProcesses(i))]),'Value')
                    userfcn_lampSwitch(childProcesses(i),1,handles)
               % ------------------------------------------ %
               else
                    % Turn on the childProcesses checkbox
                    userfcn_enable (childProcesses(i),'on',handles);
               end
            end
        % Checkbox is unselected
        case 0
            % If success = 1, release checkbox dependency enable/disable control
            if ~isempty(userData.crtPackage.processes_{index}) ...
                   && userData.crtPackage.processes_{index}.success_
                return;
            else
                for i =1:length(childProcesses)
                    % Check the childProcess requires the unselected
                    % process
                    if M(childProcesses(i),index)==2, continue; end
                    % Turn off and uncheck the follower checkboxes
                    userfcn_enable(childProcesses(i),'off',handles,true);
                
                    userfcn_lampSwitch(childProcesses(i),0,handles);
                end
            end
        otherwise
            error(['User-defined error: unexpected value of ''value'' property',...
                            'in checkbox object']);
     end
            
end




