function varargout = imKymoAnalysis(varargin)
%imKymoAnalysis : This is the GUI interface to 'imKymograph'. 
%                 See help imKymograph.

% Fields of 'handles' that contain user input data through the gui.
%    numImages     : The number of images to be analyzed.
%    imgFileList   : A cell array of the list of image files to be analyzed.
%    numKymoCurves : The number of curves or lines where the kymograph analysis
%                    is performed.
%    x (or y)      : A cell array each element of which is a numerical array
%                    that defines the x (or y) coordinates of a list of
%                    interpolation points on one curve or line. The length of x
%                    and y must be the same. And, if the length is 2, it defines
%                    the line between the two points.
%    kCurveMP      : The point we use to mark each kymo curve. It is the 
%                    left-most point of the curve.
%    width         : A cell array each element of which is the width of the
%                    tube around each curve (or line) for kymograph analysis.
%                    It has to be an odd number. If the input is even, it is
%                    added by 1.
%    selKymoCurve  : The index of the selected curve whose kymograph is to be
%                    shown and analyzed.
%    kymo          : A cell array of length 'numKymoCurves'. It stores the 
%                    kymograph images.
%    vLineX 
%    vLineY        : Coordinates of the lines drawn on the shown kymograph
%                    to estimate the velocity of flow.
%    numVLines     : The number of lines drawn on the shown kymograph
%                    for the estimation of velocity.
%    selVLine      : The selected velocity line on the kymograph.
%    whatIsShown   : A string that indicates what is shown in the figure
%                    window. Possible values: 'image' (default) and 'kymo'.
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
handles.numImages     = 0;
handles.imgFileList   = {};
handles.numKymoCurves = 0;
handles.kymoX         = {};
handles.kymoY         = {};
handles.kCurveMP      = [];
handles.width         = [];
handles.selKymoCurve  = 0;
handles.defaultWidth  = 5;
handles.whatIsShown   = 'image';
handles.numVLines     = 0;
handles.kymo          = {};
handles.vLineX        = {};
handles.vLineY        = {};
handles.vLineMP       = [];
handles.selVLine      = 0;
handles.manV          = [];
handles.calV          = [];

%Get the handles to some GUI objects.
handles.defWidthH  = findobj('tag','defaultWidth');
handles.manVFieldH = findobj('tag','manVField');
handles.calVFieldH = findobj('tag','calVField');
set(handles.defWidthH,'string',num2str(handles.defaultWidth));

%Caculate axes postion
axesWin = findobj('tag','axesWin');
handles.imgAP  = get(axesWin,'Position');
%delete(axesWin);

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

imgFileList = getFileStackNames([pathName firstImgFileName]);

%Display the first image
image = imread(imgFileList{1});
imshow(image,[]);
set(gca,'Units','pixels','Position',handles.imgAP);

numImages   = inputdlg('Number of images:','Enter Number of Images', ...
   1,{num2str(length(imgFileList))});
numImages   = str2num(numImages{1});

handles.numImages   = numImages;
handles.imgFileList = imgFileList(1:numImages);
handles.image       = image;

%Whenever new images are selected, clear the previous selected curves.
handles.numKymoCurves = 0;
handles.kymoX         = {};
handles.kymoY         = {};
handles.kCurveMP      = [];
handles.width         = [];
handles.selKymoCurve  = 0;
handles.whatIsShown   = 'image';
handles.numVLines     = 0;
handles.kymo          = {};
handles.vLineX        = {};
handles.vLineY        = {};
handles.vLineMP       = [];
handles.selVLine      = 0;
handles.manV          = [];
handles.calV          = [];

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

numKymoCurves = handles.numKymoCurves+1;

x = str2num(ans{1});
y = str2num(ans{2});
w = str2num(ans{3});

if mod(w,2) == 0
   w = w+1;
end

%Find the mark point which is the left-most end of the line.
if x(end) > x(1)
   kCurveMP = [x(1) y(1)];
else
   kCurveMP = [x(end) y(end)];
end

%Plot the curve or line.
%drawCurve(x,y,w,kCurveMP,num2str(numKymoCurves));

%Update GUI data.
handles.kymoX{numKymoCurves}      = x;
handles.kymoY{numKymoCurves}      = y;
handles.width(numKymoCurves)      = w;
handles.kCurveMP(numKymoCurves,:) = kCurveMP;
handles.numKymoCurves             = numKymoCurves;
handles.selKymoCurve              = numKymoCurves;

handles.numVLines(numKymoCurves) = 0;
handles.kymo{numKymoCurves}      = {};
handles.kymAP(numKymoCurves,:)   = zeros(1,4);
handles.vLineX{numKymoCurves}    = {};
handles.vLineY{numKymoCurves}    = {};
handles.vLineMP{numKymoCurves}   = [];
handles.selVLine(numKymoCurves)  = 0;
handles.manV{numKymoCurves}      = [];
handles.calV(numKymoCurves)      = NaN;

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

[x,y] = imSelCurve;

numKymoCurves = handles.numKymoCurves+1;

ans = inputdlg({'Width of the new line:'}, ...
   'Enter the Width',1,{num2str(handles.defaultWidth)});
w = str2num(ans{1});
if mod(w,2) == 0
   w = w+1;
end

%Find the mark point which is the left-most end of the line.
if x(end) > x(1)
   kCurveMP = [x(1) y(1)];
else
   kCurveMP = [x(end) y(end)];
end

%Plot the curve or line.
%drawCurve(x,y,w,kCurveMP,num2str(numKymoCurves));

%Update GUI data.
handles.kymoX{numKymoCurves}      = x;
handles.kymoY{numKymoCurves}      = y;
handles.width(numKymoCurves)      = w;
handles.kCurveMP(numKymoCurves,:) = kCurveMP;
handles.numKymoCurves             = numKymoCurves;
handles.selKymoCurve              = numKymoCurves;

handles.numVLines(numKymoCurves) = 0;
handles.kymAP(numKymoCurves,:)   = zeros(1,4);
handles.vLineX{numKymoCurves}    = {};
handles.vLineY{numKymoCurves}    = {};
handles.vLineMP{numKymoCurves}   = [];
handles.selVLine(numKymoCurves)  = 0;
handles.manV{numKymoCurves}      = [];
handles.calV(numKymoCurves)      = NaN;

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

numKymoCurves = handles.numKymoCurves;
if numKymoCurves == 0
   return;
end

showKymograph_Callback(hObject,eventdata,handles);
handles = guidata(hObject);

[x,y] = imSelCurve(2);

selKymoCurve = handles.selKymoCurve;
numVLines = handles.numVLines(selKymoCurve)+1;
%Find the mark point which is the left-most end of the line.
if x(end) > x(1)
   vLineMP = [x(1) y(1)];
else
   vLineMP = [x(end) y(end)];
end

%Plot the curve or line.
%drawCurve(x,y,w,kCurveMP,num2str(numKymoCurves));

%Calculate the velocity from the drawn line slope.
if abs(y(end)-y(1)) < 1
   handles.manV{selKymoCurve}(numVLines) = Inf;
else
   handles.manV{selKymoCurve}(numVLines) = (x(end)-x(1))/(y(end)-y(1))* ...
      handles.width(handles.selKymoCurve);
end

%Update GUI data.
handles.vLineX{selKymoCurve}{numVLines}    = x;
handles.vLineY{selKymoCurve}{numVLines}    = y;
handles.vLineMP{selKymoCurve}(numVLines,:) = vLineMP;
handles.numVLines(selKymoCurve)            = numVLines;
handles.selVLine(selKymoCurve)             = numVLines;

redrawAllKymo(handles);
guidata(hObject,handles);

% --- Executes on button press of 'calVelocity'.
function calVelocity_Callback(hObject, eventdata, handles)
% hObject    handle to calVelocity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

numKymoCurves = handles.numKymoCurves;
if numKymoCurves == 0
   return;
end

showKymograph_Callback(hObject,eventdata,handles);
handles = guidata(hObject);

selKymoCurve = handles.selKymoCurve;
kymo         = handles.kymo{selKymoCurve};
bw           = handles.width(selKymoCurve);

handles.calV(selKymoCurve) = imKymoSpeed(kymo,bw);
%if isnan(handles.calV(selKymoCurve))
%   handles.calV(selKymoCurve) = imKymoSpeed(kymo,bw);
%end

set(handles.calVFieldH,'String',num2str(handles.calV(selKymoCurve)));
guidata(hObject,handles);

% --- Executes on button press in showKymograph.
function showKymograph_Callback(hObject, eventdata, handles)
% hObject    handle to showKymograph (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.numKymoCurves == 0
   return;
end

if strcmp(handles.whatIsShown,'kymo') == 1
   return;
end

selKymoCurve = handles.selKymoCurve;
if selKymoCurve > length(handles.kymo) | ...
   isempty(handles.kymo{selKymoCurve})
   hmsg = msgbox('Kymograph stacking in progress ...');
   drawnow;
   handles.kymo{selKymoCurve} = imKymograph(handles.imgFileList, ...
   handles.kymoX{selKymoCurve},handles.kymoY{selKymoCurve}, ...
   handles.width(selKymoCurve));

   %Calculate the position of the axis for kymograph.
   imgAP = handles.imgAP;
   imgW  = imgAP(3); %Width of the image.
   imgH  = imgAP(4); %Height of the image.
   kymSz = size(handles.kymo{selKymoCurve});
   kymW  = min(imgW,kymSz(2)); %Width of the kymograph image.
   kymH  = min(imgH,kymSz(1)); %Height of the kymograph image.

   handles.kymAP(selKymoCurve,:) = [imgAP(1:2)+[(imgW-kymW)/2 ...
   (imgH-kymH)/2] kymW kymH];

   close(hmsg);
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
   numLines = handles.numKymoCurves;
   lineMP   = handles.kCurveMP;
elseif strcmp(handles.whatIsShown,'kymo') == 1
   selKymoCurve = handles.selKymoCurve;
   numLines = handles.numVLines(selKymoCurve);
   lineMP   = handles.vLineMP{selKymoCurve};
else
   set(gcf,'Pointer','arrow');
   return;
end

if numLines == 0
   set(gcf,'Pointer','arrow');
   return;
end

%Get the current pointer location.
p = get(gca,'CurrentPoint');

%Calculate the distance between the current pointer 'p' to the left-most end
% of each line. If the minimum distance is less than 5 pixels and the mouse
% action is left click, choose the corresponding line by highlight.
dist = sqrt((p(1,1)-lineMP(:,1)).^2+(p(1,2)-lineMP(:,2)).^2);

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

whatIsShown = handles.whatIsShown;

if strcmp(whatIsShown,'image') == 1
   selLine  = handles.selKymoCurve;
   numLines = handles.numKymoCurves;
   lineMP   = handles.kCurveMP;
elseif strcmp(whatIsShown,'kymo') == 1
   selKymoCurve = handles.selKymoCurve;
   selLine  = handles.selVLine(selKymoCurve);
   numLines = handles.numVLines(selKymoCurve);
   lineMP  = handles.vLineMP{selKymoCurve};
else
   return;
end

if numLines == 0
   return;
end

%Get the current pointer location.
p = get(gca,'CurrentPoint');

%Calculate the distance between the current pointer 'p' to the left-most end
% of each line. If the minimum distance is less than 5 pixels and the mouse
% action is left click, choose the corresponding line by highlight.
dist = sqrt((p(1,1)-lineMP(:,1)).^2+(p(1,2)-lineMP(:,2)).^2);

[minD,index] = min(dist);
if minD <= 10
   mouseAction = get(gcf,'SelectionType');
   if strcmp(mouseAction,'normal') == 1 | strcmp(mouseAction,'open') == 1
      %We have a left-mouse button click or double click. 
      % Select this line and hightlight it.
      if strcmp(whatIsShown,'image') == 1
         if index ~= selLine
            handles.selKymoCurve = index;
            redrawAllImg(handles);
         end

         if strcmp(mouseAction,'open') == 1
            %Open a dialog to edit the line.
         end
      elseif strcmp(whatIsShown,'kymo') == 1
         if index ~= selLine
            handles.selVLine(selKymoCurve) = index;
            redrawAllKymo(handles);
         end
      end
   elseif strcmp(mouseAction,'alt') == 1
      %We have a right click. Delete the corresponding line in focus and
      % redraw everything.
      if strcmp(whatIsShown,'image') == 1
         handles.numKymoCurves     = numLines-1;
         handles.kymoX(index)      = [];
         handles.kymoY(index)      = [];
         handles.width(index)      = [];
         handles.kCurveMP(index,:) = [];

         if index <= length(handles.kymo)
            handles.kymo(index) = [];
         end

         handles.numVLines(index) = [];
         handles.vLineX(index)    = [];
         handles.vLineY(index)    = [];
         handles.vLineMP(index)   = [];
         handles.selVLine(index)  = [];
         handles.manV(index)      = [];
         handles.calV(index)      = [];

         if numLines == 1
            handles.selKymoCurve = 0;
         elseif index < selLine
            handles.selKymoCurve = selLine-1;
         elseif selLine == numLines
            handles.selKymoCurve = handles.numKymoCurves;
         end

         redrawAllImg(handles);
      elseif strcmp(handles.whatIsShown,'kymo') == 1
         handles.numVLines(selKymoCurve)       = numLines-1;
         handles.vLineX{selKymoCurve}(index)   = [];
         handles.vLineY{selKymoCurve}(index)   = [];
         handles.vLineMP{selKymoCurve}(index,:) = [];
         handles.manV{selKymoCurve}(index)     = [];

         if numLines == 1
            handles.selVLine(selKymoCurve) = 0;
         elseif index < selLine
            handles.selVLine(selKymoCurve) = selLine-1;
         elseif selLine == numLines
            handles.selVLine(selKymoCurve) = numLines-1;
         end

         redrawAllKymo(handles);
      end
   end
end

guidata(hObject,handles);

% --- Some user defined subfunctions.
function drawCurve(x,y,w,kCurveMP,label)
%This function draws a curve though points defined in (x,y) and draw two bars
% that are perpendicular to the curve at the beginning and end points with
% width '2*w'. It also put a label at 'kCurveMP'. The function returns a column
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

%Label the number of the line at 'kCurveMP'. To figure out an appropriate
% distance to label, we need to scale the x and y range of the image to units
% in pixels.
set(gca,'Units','pixels');
xlim = get(gca,'XLim');
ylim = get(gca,'YLim');
axesP = get(gca,'Position');
axesW = xlim(2)-xlim(1); %The width of the axes in plotting units.
axesH = ylim(2)-ylim(1);
axesWPix = axesP(3); %The width of the axes in pixels.
axesHPix = axesP(4); %The height of the axes in pixels.
tH = text(kCurveMP(1)-10*axesW/axesWPix,kCurveMP(2)-5*axesH/axesHPix,label);
set(tH,'color','g');


function redrawAllImg(handles)
%Redraw the image and all the chosen lines (or curves).

imgAP = handles.imgAP;
delete(gca);
imshow(handles.image,[]); 
set(gca,'Units','pixels','Position',imgAP);
axis on; hold on;

for k = 1:handles.numKymoCurves
   drawCurve(handles.kymoX{k},handles.kymoY{k},floor(handles.width(k)/2), ...
      handles.kCurveMP(k,:),num2str(k));
end

index = handles.selKymoCurve;
if index ~= 0
   h = plot(handles.kymoX{index},handles.kymoY{index},'r');
   set(h,'LineWidth',1);
end

function drawKymoLine(x,y,lineMP,label)
%Draw a line on the kymograph image so that the velocity can be estimated.

h = line(x,y);
set(h,'color','g');

plot(x,y,'b.');
%Label the number of the line at 'lineMP'. To figure out an appropriate
% distance to label, we need to scale the x and y range of the image to units
% in pixels.
set(gca,'Units','pixels');
xlim = get(gca,'XLim');
ylim = get(gca,'YLim');
axesP = get(gca,'Position');
axesW = xlim(2)-xlim(1); %The width of the axes in plotting units.
axesH = ylim(2)-ylim(1);
axesWPix = axesP(3); %The width of the axes in pixels.
axesHPix = axesP(4); %The height of the axes in pixels.
tH = text(lineMP(1)-10*axesW/axesWPix,lineMP(2)-5*axesH/axesHPix,label);
set(tH,'color','g');


function redrawAllKymo(handles)
%Redraw the kymograph for the selected curve and any manual tracking velocity
% line.

selKymoCurve = handles.selKymoCurve;

%handles.kymAP stores the axes position of the kymograph.
kymAP = handles.kymAP(selKymoCurve,:);

delete(gca);
axes('Units','pixels','Position',kymAP);
imshow(handles.kymo{selKymoCurve},[]);
hold on;
bandW     = handles.width(selKymoCurve);
numImages = handles.numImages;
if numImages > 3
   ytickL = [1 ceil(numImages/2) numImages];
else
   ytickL = [1 numImages];
end
ytick = bandW*ytickL;
set(gca,'YTick',ytick,'YTickLabel',ytickL);
axis on;
title(['Kymograph ' num2str(selKymoCurve)]);
xlabel('pixel');
ylabel('frame');

numVLines = handles.numVLines(selKymoCurve);
if numVLines == 0
   set(handles.manVFieldH,'String','');

   if isnan(handles.calV(selKymoCurve))
      set(handles.calVFieldH,'String','');
   else
      set(handles.calVFieldH,'String',num2str(handles.calV(selKymoCurve)));
   end
   
   return;
end

for j = 1:numVLines
   drawKymoLine(handles.vLineX{selKymoCurve}{j}, ...
      handles.vLineY{selKymoCurve}{j}, handles.vLineMP{selKymoCurve}(j,:), ...
      num2str(j));
end

index = handles.selVLine(selKymoCurve);
if index ~= 0
   h = line(handles.vLineX{selKymoCurve}{index}, ...
      handles.vLineY{selKymoCurve}{index});
   set(h,'color','r','LineWidth',1);

   set(handles.manVFieldH,'String',num2str(handles.manV{selKymoCurve}(index)));
else
   set(handles.manVFieldH,'String','');
end

if isnan(handles.calV(selKymoCurve))
   set(handles.calVFieldH,'String','');
else
   set(handles.calVFieldH,'String',num2str(handles.calV(selKymoCurve)));
end
