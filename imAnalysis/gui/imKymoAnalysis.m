function varargout = imKymoAnalysis(varargin)
%imKymoAnalysis : This is the GUI interface to 'imKymograph'. 
%                 See help imKymograph.
%
% Fields of 'handles' that contain user input data through the gui.
%    numImages   : The number of images to be analyzed.
%    imgFileList : A cell array of the list of image files to be analyzed.
%    numCurves    : The number of curves or lines where the kymograph analysis
%                  is performed.
%    x (or y)    : A cell array each element of which is a numerical array
%                  that defines the x (or y) coordinates of a list of
%                  interpolation points on one curve or line. The length of x
%                  and y must be the same. And, if the length is 2, it defines
%                  the line between the two points.
%    markP       : The point we use to mark each line. It is the left-most
%                  point of the line.
%    width       : A cell array each element of which is the width of the
%                  section around one curve (or line).
%
% IMKYMOANALYSIS M-file for imKymoAnalysis.fig
%      IMKYMOANALYSIS, by itself, creates a new IMKYMOANALYSIS or raises the existing
%      singleton*.
%
%      H = IMKYMOANALYSIS returns the handle to a new IMKYMOANALYSIS or the handle to
%      the existing singleton*.
%
%      IMKYMOANALYSIS('Property','Value',...) creates a new IMKYMOANALYSIS using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to imKymoAnalysis_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      IMKYMOANALYSIS('CALLBACK') and IMKYMOANALYSIS('CALLBACK',hObject,...) call the
%      local function named CALLBACK in IMKYMOANALYSIS.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help imKymoAnalysis

% Last Modified by GUIDE v2.5 02-Apr-2004 15:59:23

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @imKymoAnalysis_OpeningFcn, ...
                   'gui_OutputFcn',  @imKymoAnalysis_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
                   'gui_Callback',   []);
if nargin & isstr(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before imKymoAnalysis is made visible.
function imKymoAnalysis_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for imKymoAnalysis
handles.output = hObject;

%Initialization
handles.numImages    = 0;
handles.imgFileList  = {};
handles.numCurves    = 0;
handles.x            = {};
handles.y            = {};
handles.markP        = [];
handles.width        = {};
handles.selCurve     = 0;
handles.defaultWidth = 10;
handles.whatIsShown  = 'none';
handles.numKymoLines = 0;
handles.kymo         = {};
handles.kymoX        = {};
handles.kymoY        = {};
handles.kymoMP       = [];
handles.selKymoLine  = 0;

h = findobj('tag','defaultWidth');
set(h,'string',num2str(handles.defaultWidth));

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes imKymoAnalysis wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = imKymoAnalysis_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in selectImages.
function selectImages_Callback(hObject, eventdata, handles)
% hObject    handle to selectImages (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[firstImgFileName pathName filterIndex] = uigetfile({'*.tif;*.gif;*.jpg', ...
   'Image Files (*.tif,*.gif,*.jpg)';
   '*.tif','TIFF files (*.tif)';
   '*.gif','GIF files (*.gif)';
   '*.jpg','JPEG files (*.jpg)'},'Select the first image');

if filterIndex == 0
   return;
end

imgFileList = getFileStackNames(firstImgFileName);

%Display the first image
image = imread(imgFileList{1});
cla;
imshow(image,[]);

numImages   = inputdlg('Number of images:','Enter Number of Images', ...
   1,{num2str(length(imgFileList))});
numImages   = str2num(numImages{1});

handles.numImages   = numImages;
handles.imgFileList = imgFileList(1:numImages);
handles.image       = image;

%Whenever new images are selected, clear the previous selected curves.
handles.numCurves    = 0;
handles.x            = {};
handles.y            = {};
handles.markP        = [];
handles.width        = {};
handles.selCurve     = 0;
handles.whatIsShown  = 'image';
handles.numKymoLines = 0;
handles.kymo         = {};
handles.kymoX        = {};
handles.kymoY        = {};
handles.kymoMP       = [];
handles.selKymoLine  = 0;

%Display the first image
redrawAllImg(handles);
guidata(hObject,handles);

% --- Executes on button press in lineDialog.
function lineDialog_Callback(hObject, eventdata, handles)
% hObject    handle to lineDialog (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if strcmp(handles.whatIsShown,'image') == 0 | handles.numImages == 0
   %When what is shown in the figure window is not the image or when no image
   % is input yet.
   return;
end

ans = inputdlg({'X coordinates:','Y coordinates:','Width:'}, ...
   'Enter a New Line',1,{'' '' num2str(handles.defaultWidth)});

numCurves = handles.numCurves+1;

x = str2num(ans{1});
y = str2num(ans{2});
w = str2num(ans{3});

%Find the mark point which is the left-most end of the line.
if x(end) > x(1)
   markP = [x(1) y(1)];
else
   markP = [x(end) y(end)];
end

%Plot the curve or line.
%drawCurve(x,y,w,markP,num2str(numCurves));

%Update GUI data.
handles.x{numCurves}       = x;
handles.y{numCurves}       = y;
handles.width{numCurves}   = w;
handles.markP(numCurves,:) = markP;
handles.numCurves          = numCurves;
handles.selCurve           = numCurves;

handles.numKymoLines = 0;
handles.kymo         = {};
handles.kymoX        = {};
handles.kymoY        = {};
handles.kymoMP       = [];
handles.selKymoLine  = 0;

redrawAllImg(handles);
guidata(hObject,handles);

% --- Executes on button press in lineDraw.
function lineDraw_Callback(hObject, eventdata, handles)
% hObject    handle to lineDraw (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if strcmp(handles.whatIsShown,'image') == 0 | handles.numImages == 0
   %When what is shown in the figure window is not the image or when no image
   % is input yet.
   return;
end

[x,y] = roicurve;

numCurves = handles.numCurves+1;

ans = inputdlg({'Width of the new line:'}, ...
   'Enter the Width',1,{num2str(handles.defaultWidth)});
w = str2num(ans{1});

%Find the mark point which is the left-most end of the line.
if x(end) > x(1)
   markP = [x(1) y(1)];
else
   markP = [x(end) y(end)];
end

%Plot the curve or line.
%drawCurve(x,y,w,markP,num2str(numCurves));

%Update GUI data.
handles.x{numCurves}       = x;
handles.y{numCurves}       = y;
handles.width{numCurves}   = w;
handles.markP(numCurves,:) = markP;
handles.numCurves          = numCurves;
handles.selCurve           = numCurves;

handles.numKymoLines = 0;
handles.kymo         = {};
handles.kymoX        = {};
handles.kymoY        = {};
handles.kymoMP       = [];
handles.selKymoLine  = 0;

redrawAllImg(handles);
guidata(hObject,handles);

% --- Executes on the text field 'defaultWidth'.
function defaultWidth_Callback(hObject, eventdata, handles)
% hObject    handle to defaultWidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
entry = get(hObject,'string');
if isnan(str2double(entry))
   errordlg('You must enter a numeric value','Bad Input','modal');
end

handles.defaultWidth = str2double(entry);
guidata(hObject,handles);

% --- Executes on button press in showImage.
function showImage_Callback(hObject, eventdata, handles)
% hObject    handle to showImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if strcmp(handles.whatIsShown,'image') == 1
   return;
end

redrawAllImg(handles);
handles.whatIsShown = 'image';

guidata(hObject,handles);

% --- Executes on button press of 'manVelocityTrack'.
function manVelocityTrack_Callback(hObject, eventdata, handles)
% hObject    handle to manVelocityTrack (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if strcmp(handles.whatIsShown,'kymo') == 0
   %When what is shown in the figure window is not the kymograph.
   return;
end

[x,y] = roicurve;

numKymoLines = handles.numKymoLines+1;

%Find the mark point which is the left-most end of the line.
if x(end) > x(1)
   markP = [x(1) y(1)];
else
   markP = [x(end) y(end)];
end

%Plot the curve or line.
%drawCurve(x,y,w,markP,num2str(numCurves));

%Update GUI data.
handles.kymoX{numKymoLines}    = x;
handles.kymoY{numKymoLines}    = y;
handles.kymoMP(numKymoLines,:) = markP;
handles.numKymoLines           = numKymoLines;
handles.selKymoLine            = numKymoLines;

redrawAllKymo(handles);
guidata(hObject,handles);


% --- Executes on button press in showKymograph.
function showKymograph_Callback(hObject, eventdata, handles)
% hObject    handle to showKymograph (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.numKymoLines == 0
   return;
end

if strcmp(handles.whatIsShown,'kymo') == 1
   return;
end

if handles.selCurve > length(handles.kymo) | ...
   isempty(handles.kymo{handles.selCurve})
   handles.kymo{selCurve} = imKymograph(handles.imgFileList, ...
      handles.x{selCurve},handles.y{selCurve}, ...
      handles.width{selCurve});
end

redrawAllKymo(handles);
handles.whatIsShown = 'kymo';

guidata(hObject,handles);


% --- Executes on mouse button motion.
function winButtonMotion_Callback(hObject, eventdata, handles)
% hObject    handle to 'WindowButtonMotionFcn' (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if strcmp(handles.whatIsShown,'image') == 1
   numCurves = handles.numCurves;
   markP     = handles.markP;
elseif strcmp(handles.whatIsShown,'kymo') == 1
   numCurves = handles.numKymoLines;
   markP     = handles.kymoMP;
else
   set(gcf,'Pointer','arrow');
   return;
end

if numCurves == 0
   set(gcf,'Pointer','arrow');
   return;
end

%Get the current pointer location.
p = get(gca,'CurrentPoint');

%Calculate the distance between the current pointer 'p' to the left-most end
% of each line. If the minimum distance is less than 5 pixels and the mouse
% action is left click, choose the corresponding line by highlight.
dist = sqrt((p(1,1)-markP(:,1)).^2+(p(1,2)-markP(:,2)).^2);

[minD,index] = min(dist);
if minD <= 10
   %Change the mouse shape.
   set(gcf,'Pointer','circle');
else
   set(gcf,'Pointer','arrow');
end

% --- Executes on mouse button down.
function winButtonDown_Callback(hObject, eventdata, handles)
% hObject    handle to 'WindowButtonDownFcn' (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if strcmp(handles.whatIsShown,'image') == 1
   numCurves = handles.numCurves;
   markP     = handles.markP;
elseif strcmp(handles.whatIsShown,'kymo') == 1
   numCurves = handles.numKymoLines;
   markP     = handles.kymoMP;
else
   return;
end

if numCurves == 0
   return;
end

%Get the current pointer location.
p = get(gca,'CurrentPoint');

%Calculate the distance between the current pointer 'p' to the left-most end
% of each line. If the minimum distance is less than 5 pixels and the mouse
% action is left click, choose the corresponding line by highlight.
dist = sqrt((p(1,1)-markP(:,1)).^2+(p(1,2)-markP(:,2)).^2);

[minD,index] = min(dist);
if minD <= 10
   mouseAction = get(gcf,'SelectionType');
   if strcmp(mouseAction,'normal') == 1 | strcmp(mouseAction,'open') == 1
      %We have a left-mouse button click or double click. 
      % Select this line and hightlight it.
      if strcmp(handles.whatIsShown,'image') == 1
         handles.selCurve = index;
         cla; redrawAllImg(handles);
         if strcmp(mouseAction,'open') == 1
            %Open a dialog to edit the line.
         end
      elseif strcmp(handles.whatIsShown,'kymo') == 1
         handles.selKymoLine = index;
         cla; redrawAllKymo(handles);
      end
   elseif strcmp(mouseAction,'alt') == 1
      %We have a right click. Delete the corresponding line in focus and
      % redraw everything.
      if strcmp(handles.whatIsShown,'image') == 1
         handles.numCurves      = handles.numCurves-1;
         handles.x(index)       = [];
         handles.y(index)       = [];
         handles.width(index)   = [];
         handles.markP(index,:) = [];
         handles.selCurve       = handles.numCurves;
         handles.kymo{index}    = [];
         redrawAllImg(handles);
      elseif strcmp(handles.whatIsShown,'kymo') == 1
         handles.numKymoLines    = handles.numCurves-1;
         handles.kymoX(index)    = [];
         handles.kymoY(index)    = [];
         handles.width(index)    = [];
         handles.kymoMP(index,:) = [];
         handles.selKymoLine     = handles.numKymoLines;
         redrawAllKymo(handles);
      end
   end
end

guidata(hObject,handles);

% --- Some user defined subfunctions.
function drawCurve(x,y,w,markP,label)
%This function draws a curve though points defined in (x,y) and draw two bars
% that are perpendicular to the curve at the beginning and end points with
% width '2*w'. It also put a label at 'markP'. The function returns a column
% vector of handles to all the drawn objects.

plot(x,y,'g');

%Draw a bar perpendicular to the curve with width '2*w' at the first point.
len = sqrt((x(2)-x(1))^2+(y(2)-y(1))^2);
xL  = x(1) - w*(y(1)-y(2))/len;
yL  = y(1) - w*(x(2)-x(1))/len;
xR  = x(1) + w*(y(1)-y(2))/len;
yR  = y(1) + w*(x(2)-x(1))/len;
plot([xL xR],[yL yR],'y');

%Draw a bar perpendicular to the curve with width '2*w' at the last point.
len = sqrt((x(end)-x(end-1))^2+(y(end)-y(end-1))^2);
xL  = x(end) - w*(y(end-1)-y(end))/len;
yL  = y(end) - w*(x(end)-x(end-1))/len;
xR  = x(end) + w*(y(end-1)-y(end))/len;
yR  = y(end) + w*(x(end)-x(end-1))/len;
plot([xL xR],[yL yR],'y');
plot(x,y,'bo');

%Label the number of the new curve.
tH = text(markP(1)-25,markP(2)-5,label);
set(tH,'color','g');


function redrawAllImg(handles)
%Redraw the image and all the chosen lines (or curves).

cla;
imshow(handles.image,[]); axis on; hold on;

for k = 1:handles.numCurves
   drawCurve(handles.x{k},handles.y{k},handles.width{k}, ...
      handles.markP(k,:),num2str(k));
end

index = handles.selCurve;
if index ~= 0
   h = plot(handles.x{index},handles.y{index},'r');
   set(h,'LineWidth',2);
end

function drawKymoLine(x,y,markP,label)
%Draw a line on the kymograph image so that the velocity can be estimated.

h = line(x,y);
set(h,'color','g');

plot(x,y,'b.');
%Label the number of the line.
tH = text(markP(1)-25,markP(2)-5,label);
set(tH,'color','g');


function redrawAllKymo(handles)
%Redraw the kymograph for the selected curve and any manual tracking velocity
% line.

cla;
imshow(handles.kymo{selCurve},[]); hold on;

numKymoLines = handles.numKymoLines;
if numKymoLines == 0
   return;
end

for k = 1:numKymoLines
   drawKymoLine(handles.kymoX{k},handles.kymoY{k},handles.kymoMP(k,:), ...
      num2str(k));
end

index = handles.selKymoLine;
if index ~= 0
   h = line(handles.kymoX{index},handles.kymoY{index});
   set(h,'color','r','LineWidth',1);
end
