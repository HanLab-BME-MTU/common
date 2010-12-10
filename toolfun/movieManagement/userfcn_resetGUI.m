function userfcn_resetGUI(handles)
% GUI tool function: reset the GUI
% 
% Input: 
%       handles - the "handles" of package GUI control panel
%
%
% Chuangang Ren
% 08/2010

userData = get(handles.figure1, 'UserData');
l = size(userData.dependM,1);
userfcn_drawIcon(handles, 'clear', 1:l);
userfcn_enable(1:l, 'on', handles)

for i = 1:l
    eval(['set(handles.checkbox_' num2str(i) ', ''FontWeight'', ''normal'', ''Value'', 0) '])
    eval(['set(handles.pushbutton_show_' num2str(i) ', ''Enable'', ''off'') '])
end