function [xi,yi] = roicurve(varargin)
%roicurve: Interactively insert a curve in the current figure without 
%          acturally drawing it.
%
% SYNTAX: [xi,yi]= roicurve
%    The coordinates of the curve is returned in 'xi' and 'yi'.

if nargin > 0
   error('No input argument is implemented yet.');
end

%Get the current figure.
h = gcf; figure(h);

%Initialize with empty data.
handles.numPoints = 0;
handles.xi        = [];
handles.yi        = [];
handles.dashL     = {};
handles.dotL      = {};
handles.mDashL    = [];
handles.mDotL     = [];

%Update GUI data.
guidata(h,handles);

%Get old properties.
oldPointer = get(h,'Pointer');
oldButtonDownFcn = get(h,'WindowButtonDownFcn');
oldButtonMotionFcn = get(h,'WindowButtonMotionFcn');
oldKeyPressFcn = get(h,'KeyPressFcn');

set(h,'DoubleBuffer','on','Pointer','cross');
set(h,'WindowButtonDownFcn',@buttonDown_Callback);
set(h,'WindowButtonMotionFcn',@buttonMove_Callback);
set(h,'KeyPressFcn',@keyPress_Callback);

%Set it to be modal.
%set(h,'WindowStyle','modal','MenuBar','figure');

w    = waitforbuttonpress;
ckey = get(h,'CurrentCharacter');
while w ~= 1 | (w == 1 & double(ckey) ~= 13)
   %When there is no key stroke from the key board or the key stroke is not
   % return, keep waiting.
   w = waitforbuttonpress;
   ckey = get(h,'CurrentCharacter');
end

%Set it back to normal.
set(h,'WindowStyle','normal','Pointer',oldPointer);
set(h,'WindowButtonDownFcn',oldButtonDownFcn);
set(h,'WindowButtonMotionFcn',oldButtonMotionFcn);
set(h,'KeyPressFcn',oldKeyPressFcn);

handles = guidata(h);

xi = handles.xi;
yi = handles.yi;

%Delete the lines that are drawn during selection.
if handles.numPoints > 1
   for k = 1:handles.numPoints-1
      delete(handles.dashL{k});
      delete(handles.dotL{k});
   end
end

if ~isempty(handles.mDashL)
   delete(handles.mDashL);
   delete(handles.mDotL);
end


% --- Call back function for WindowButtonDownFcn.
function buttonDown_Callback(obj,eventdata)

%Get the GUI data
handles = guidata(obj);

numPoints = handles.numPoints+1;

%Delete previously drawn line while mounse is moving.
if ~isempty(handles.mDashL)
   delete(handles.mDashL);
   delete(handles.mDotL);
   handles.mDashL = [];
   handles.mDotL = [];
end

p = get(gca,'CurrentPoint');

if numPoints > 1
   handles.dashL{numPoints-1} = line([handles.xi(end);p(1,1)], ...
      [handles.yi(end);p(1,2)]);
   set(handles.dashL{numPoints-1},'LineStyle','--');
   set(handles.dashL{numPoints-1},'color','k');
   handles.dotL{numPoints-1} = line([handles.xi(end);p(1,1)], ...
      [handles.yi(end);p(1,2)]);
   set(handles.dotL{numPoints-1},'LineStyle',':');
   set(handles.dotL{numPoints-1},'color','w');
end

handles.xi(numPoints) = p(1,1);
handles.yi(numPoints) = p(1,2);

handles.numPoints = numPoints;

%Update GUI data.
guidata(obj,handles);

% --- Call back function for WindowButtonMotionFcn.
function buttonMove_Callback(obj,eventdata)

%Get the GUI data
handles = guidata(obj);

if handles.numPoints == 0
   return;
end

%Delete previously drawn line while mounse is moving.
if ~isempty(handles.mDashL)
   delete(handles.mDashL);
   delete(handles.mDotL);
end

p = get(gca,'CurrentPoint');

%Draw the moving line
handles.mDashL = line([handles.xi(end);p(1,1)],[handles.yi(end);p(1,2)]);
set(handles.mDashL,'LineStyle','--');
set(handles.mDashL,'color','k');
handles.mDotL = line([handles.xi(end);p(1,1)],[handles.yi(end);p(1,2)]);
set(handles.mDotL,'LineStyle',':');
set(handles.mDotL,'color','w');

%Update GUI data.
guidata(obj,handles);

% --- Call back function for KeyPressFcn.
function keyPress_Callback(obj,eventdata)

handles = guidata(obj);

numPoints = handles.numPoints;
if numPoints == 0
   return;
end

p    = get(gca,'CurrentPoint');
ckey = get(gcf,'CurrentCharacter');

if double(ckey) == 8
   %When the key stroke is backspace (ascii code 8), delete the last selected 
   % point.
   handles.xi(end)   = [];
   handles.yi(end)   = [];
   handles.numPoints = numPoints-1;

   %Delete the last line.
   if numPoints > 1
      delete(handles.dashL{end});
      delete(handles.dotL{end});
      handles.dashL(end) = [];
      handles.dotL(end)  = [];
   end

   if ~isempty(handles.mDashL)
      delete(handles.mDashL);
      delete(handles.mDotL);
      handles.mDashL = [];
      handles.mDotL  = [];
   end

   if numPoints > 1
      handles.mDashL = line([handles.xi(end);p(1,1)],[handles.yi(end);p(1,2)]);
      set(handles.mDashL,'LineStyle','--');
      set(handles.mDashL,'color','k');
      handles.mDotL = line([handles.xi(end);p(1,1)],[handles.yi(end);p(1,2)]);
      set(handles.mDotL,'LineStyle',':');
      set(handles.mDotL,'color','w');
   end
end

%Update GUI data.
guidata(obj,handles);

