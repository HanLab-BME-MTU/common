function userfcn_drawIcon(handles, type, ID, msg)
% This function should be defined as a global function in future (with
% multiple packages).
%
%
userData = get(handles.figure1, 'UserData');
switch type
    case 'pass'
        iconData = userData.passIconData;
    case 'error'
        iconData = userData.errorIconData;
    case 'warn'
        iconData = userData.warnIconData;
    case 'clear'
        for i = ID
            eval( [ 'cla(handles.axes_icon' ,num2str(i), ');' ])
        end
        return
    otherwise
        error('User-defined: function input ''type'' is incorrect');
end
for i = ID
    eval(['axes(handles.axes_icon',num2str(i),');'])
    Img = image(iconData);
    set(gca, 'XLim',get(Img,'XData'),'YLim',get(Img,'YData'),...
                                        'visible','off','YDir','reverse');
    set(Img,'ButtonDownFcn',@icon_ButtonDownFcn);
    set(Img,'UserData',msg);
    eval([ 'set(Img,''tag'',''',num2str(i),''');' ])
end
