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
handles.line      = {};
handles.movingL   = {};

%Update GUI data.
guidata(h,handles);

%Get old properties.
oldPointer = get(h,'Pointer');
oldButtonDownFcn = get(h,'WindowButtonDownFcn');
oldButtonMotionFcn = get(h,'WindowButtonMotionFcn');

set(h,'DoubleBuffer','on','Pointer','cross');
set(h,'WindowButtonDownFcn',@buttonDown_Callback);
set(h,'WindowButtonMotionFcn',@buttonMove_Callback);

%Set it to be modal.
%set(h,'WindowStyle','modal','MenuBar','figure');

w = waitforbuttonpress;
while w ~= 1
   w = waitforbuttonpress;
end

%Set it back to normal.
set(h,'WindowStyle','normal','Pointer',oldPointer);
set(h,'WindowButtonDownFcn',oldButtonDownFcn);
set(h,'WindowButtonMotionFcn',oldButtonMotionFcn);

handles = guidata(h);

xi = handles.xi;
yi = handles.yi;

%Delete the lines that are drawn during selection.
if handles.numPoints > 1
   for k = 1:handles.numPoints-1
      delete(handles.line{k});
   end
end

if ~isempty(handles.movingL)
   delete(handles.movingL);
end


% --- Call back function for Window Button Down.
function buttonDown_Callback(obj,eventdata)

%Get the GUI data
handles = guidata(obj);

numPoints = handles.numPoints+1;

%Delete previously drawn line while mounse is moving.
if ~isempty(handles.movingL)
   delete(handles.movingL);
   handles.movingL = {};
end

p = get(gca,'CurrentPoint');

if numPoints > 1
   handles.line{numPoints-1} = line([handles.xi(end);p(1,1)], ...
      [handles.yi(end);p(1,2)]);
   set(handles.line{numPoints-1},'LineStyle','-.');
end

handles.xi(numPoints) = p(1,1);
handles.yi(numPoints) = p(1,2);

handles.numPoints = numPoints;

%Update GUI data.
guidata(obj,handles);

% --- Call back function for Window Button Move.
function buttonMove_Callback(obj,eventdata)

%Get the GUI data
handles = guidata(obj);

if handles.numPoints == 0
   return;
end

%Delete previously drawn line while mounse is moving.
if ~isempty(handles.movingL)
   delete(handles.movingL);
end

p = get(gca,'CurrentPoint');

%Draw the moving line
handles.movingL = line([handles.xi(end);p(1,1)],[handles.yi(end);p(1,2)]);
set(handles.movingL,'LineStyle','-.');

%Update GUI data.
guidata(obj,handles);

