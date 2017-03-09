function out = cliGUI(string,hfig1,procID,varargin)
    % Generic Command Line Interface GUI for arbitrary processes
    %
    % The user is presented with a command line interface via keyboard.
    % The command dbquit will exit debug mode without applying settings.
    % The command dbcont will apply funParams to the process being
    % modified using processGUI_ApplyFcn
    %
    % See also dbcont, dbquit, processGUI_ApplyFcn, noSettingsProcessGUI
    out = noSettingsProcessGUI(string,hfig1,procID,varargin{:});
    hObject = out;
    eventdata = [];
    handles = guihandles(out);
    set(handles.text34,'String','Check your command line interface');
    
    userData = get(out,'UserData');
    proc = userData.crtProc;
    
    % Display text
    set(handles.pushbutton_done,'Callback',@(hObject,eventdata) disp('Type <a href="matlab:dbcont">dbcont</a> to apply your settings.'));
    set(handles.pushbutton_cancel,'Callback',@(hobject,evendata) disp('Type <a href="matlab:dbquit">dbquit</a> to cancel.'));
    disp('Welcome to Command Line Interface "GUI"');
    disp('---------------------------------------');
    disp('Stored in the struct <a href="matlab:funParams">funParams</a> are the process'' parameters.');
    disp('Type <a href="matlab:dbcont">dbcont</a> to apply your changes to funParam.');
    disp('Type <a href="matlab:dbquit">dbquit</a> to cancel.')
    who
    funParams = userData.crtProc.getParameters()
    
    % Turn over command to the user
    keyboard;
    
    % Set the parameters if we get to this point via dbcont
    processGUI_ApplyFcn(hObject, eventdata, handles,funParams,varargin);
    % Delete the figure if it still exists
    if(ishandle(hObject) && isvalid(hObject))
        delete(hObject);
    end
end
