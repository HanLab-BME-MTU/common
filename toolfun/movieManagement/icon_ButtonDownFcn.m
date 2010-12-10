function icon_ButtonDownFcn(hObject, eventdata)
% This function call up a help dialog box when user click any of the icons
% in all GUIs.
%
handles = guidata(hObject);
userData = get(handles.figure1, 'UserData');
% Help dialog from MovieData panel
splitTag = regexp(get(get(hObject,'parent'), 'tag'), '_','split');

% Pass handle to userData
% If from package GUI, call pre-defined help dialog
% if called from setting GUI, call user-defined help dialog 'msgboxGUI'

if isfield(userData, 'crtProc')
    % Help dialog from setting panel
    if ~isempty(userData.crtProc)
        userData.helpFig = msgboxGUI('Text', sprintf([get(hObject,'UserData'), ...
            '\n', userfcn_copyright ]),'Title',['Help - ' userData.crtProc.name_] );
    else
        userData.helpFig = msgboxGUI('Text', sprintf([get(hObject,'UserData'), ...
            '\n', userfcn_copyright ]),'Title','Help');
    end


elseif strcmp(splitTag{1}, 'axes') && length(splitTag) >1

        if strcmpi(splitTag{2}, 'help') % Help icon
            if length(splitTag) < 3
                % Package help
                userData.packageHelpFig = msgbox(sprintf(get(hObject,'UserData')), ...
                    ['Help - ' userData.crtPackage.name_], 'custom', get(hObject,'CData'), userData.colormap, 'replace');
            else
                % Process help
                procID = str2double(splitTag{3});
                if ~isnan(procID)
                        
                    procName = regexp(userData.crtPackage.processClassNames_{procID}, 'Process','split');
                    userData.processHelpFig(procID) = msgbox(sprintf(get(hObject,'UserData')), ...
                     ['Help - ' procName{1}], 'custom', get(hObject,'CData'), userData.colormap, 'replace');
                end
            end
        
        else % Process status icon
             
            userData.iconHelpFig = msgbox(get(hObject,'UserData'), ...
                'Help', 'custom', get(hObject,'CData'), userData.colormap, 'replace');            
        end
    
else
    userData.iconHelpFig = msgbox(get(hObject,'UserData'), ...
        'Help', 'custom', get(hObject,'CData'), userData.colormap, 'replace'); 
end

set(handles.figure1, 'UserData', userData);