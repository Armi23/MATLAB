function varargout = LesionPlot(varargin)
% LESIONPLOT MATLAB code for LesionPlot.fig
%      LESIONPLOT, by itself, creates a new LESIONPLOT or raises the existing
%      singleton*.
%
%      H = LESIONPLOT returns the handle to a new LESIONPLOT or the handle to
%      the existing singleton*.
%
%      LESIONPLOT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in LESIONPLOT.M with the given input arguments.
%
%      LESIONPLOT('Property','Value',...) creates a new LESIONPLOT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before LesionPlot_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to LesionPlot_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help LesionPlot

% Last Modified by GUIDE v2.5 22-Jul-2013 11:02:15

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @LesionPlot_OpeningFcn, ...
                   'gui_OutputFcn',  @LesionPlot_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before LesionPlot is made visible.
function LesionPlot_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to LesionPlot (see VARARGIN)

handles.mode = 0; % start in length measurement mode
handles.threshold = 0.8; %for b/w boundary finding
handles.filenameindex = 1; %which of the multiple files we're on
handles.object = 1; %refers to object we're working with in slide
handles.downsample = 0.1; % default downsample amount
handles.noisecancel = 2000; % smallest object to be bounded
handles.pointmode = 'normal'; %we start in brain mode. 
handles.colorPlots = []; % Eventual output file for Length mode
handles.points = []; % Eventual output for Area mode
handles.midlinepoint = []; % Does not need counter, only allowed 2 points
handles.cpoints = []; % Points for the part of the lesion we're sure about
handles.normal = 1; % counter for drawing liberal lesion/lesion area
handles.certain = 1; % counter for when we're sure about the lesion
handles.outFile = 'default'; % default name for output file
handles.row = 1; % keep track of which row of brain we're on
handles.bregma = 0; % put bregma out of bounds first

% Choose default command line output for LesionPlot
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes LesionPlot wait for user response (see UIRESUME)
% uiwait(handles.figure1);

function d = mydistance(p1,p2)
x1 = p1(1);
x2 = p2(1);
y1 = p1(2);
y2 = p2(2);
d = sqrt((x2 - x1)^2 + (y2 - y1)^2);

function d = distancetoline(P,Q1,Q2)
d = abs(det([Q2-Q1;P-Q1]))/norm(Q2-Q1);

function pointindex = indexOfPlot(region,clickedpoint)
dist = Inf;
[m,~] = size(region);
pointindex = 0;
for i = 1:m
    check = mydistance(region(i,:),clickedpoint);
    if check < dist
        dist = check;
        pointindex = i;
    end
end

% --- Outputs from this function are returned to the command line.
function varargout = LesionPlot_OutputFcn(hObject, ~, handles)  %#ok<INUSL>
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in GenerateButton.
function GenerateButton_Callback(hObject, ~, handles) %#ok<DEFNU>
% hObject    handle to GenerateButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% File is read in and has to be made black and white to allow the 
% boundary recognition program to work. Levels of this can be changed
% in the im2bw function, but I've found that 0.8 works pretty well. Also,
% we invert it because it switches background and foreground
I = handles.image;
BW = im2bw(I,handles.threshold);
BW = ~BW; 
BW_filled = imfill(BW,'holes');
boundaries = bwboundaries(BW_filled);

% Initialized values to help scan through the boundaries we have
bsize = size(boundaries); 
i = 0; %used and incremented when making indices of real boundaries array

% Check the size of each boundary and if its small we can assume its 
% just noise or unnecessary extra stuff. We'll plot the boundaries that we
% do find to be important and attach event listeners to each point so that
% we can react if it is selected
hold on;
handles.real= [];
for k=1:bsize
    [m,~] = size(boundaries{k});
    if (m > handles.noisecancel)
                
        % columns are switched to fix weird convention of boundaries from
        % y:x to x:y. Important for plotting later. 
        i = i + 1;
        handles.real{i, 1}(:,[1,2]) = boundaries{k}(:,[2,1]);
        r = handles.real{i,1}; % real refers to our important boundaries
        
        % we keep the handle of the plots in case we want to delete later
        handles.pHandle(i) = plot(r(:,1),r(:,2),'r','LineWidth',3);
        set(handles.pHandle(i),'ButtonDownFcn',{@PointSelect,i});
        guidata(hObject,handles);
    end
end

% generate prompt for user
[rsize,~] = size(handles.real);
line = ['Generated boundaries for ' num2str(rsize) ' objects'];
PrintText(hObject,handles,line);
guidata(hObject,handles);

% --- Called when a point on a boundary is clicked. 
function PointSelect (hObject,eventdata,object)

handles = guidata(hObject); % get current handles

% In normal mode, we need to be in Area mode for boundary clicks to matter
if (strcmp(handles.pointmode, 'normal')) 
    if ~(handles.mode) % in Area mode
        
        % We find out where the mouse clicked and then use that to find the
        % index of the array that marks that point. Used distance equation. 
        coordinates = get(gca,'CurrentPoint');
        clickedpoint = coordinates(1,1:2);
        
        % Real is the important boundary that was clicked on. We find out
        % location on it by finding the closest point to our click
        % point and the index of the point becomes pointindex.
        real = handles.real{object,1};
        pointindex = indexOfPlot(real,clickedpoint);

        % We only want 2 points on the boundary. More is taken as a new set
        [~, n] = size(handles.indexset);
        if n >= 2 
            handles.indexset = [];
        end

        % Now we add the point, provide feedback about where we are.  
        handles.indexset = [handles.indexset,pointindex];
        handles.object = object; % To keep track of used bound
        set(handles.ObjectArea,'String',num2str(object));

        % Print feedback
        [~, n] = size(handles.indexset);
        line = ['Point ' num2str(n) ' selected on object' num2str(object)];
        PrintText(hObject,handles,line);
        
        if n == 2 
             
            % Find out what object we're on, get its boundaries,
            boundaryindex = handles.object;
            
            % Edit boundaries based on what points we clicked on to create
            % "boundset," and then plot this edited region. We keep track 
            % of it in case of deletion or if we want to use it more later
            boundary = handles.real{boundaryindex,1};
            boundset = cat(1,boundary(1:handles.indexset(1,1),1:2), ...
                handles.points,boundary(handles.indexset(1,2):end,1:2));
            handles.editedregion = plot(boundset(:,1),boundset(:,2),...
                'b','LineWidth',3);
            handles.boundset{boundaryindex,1} = boundset;

            % We re-initialize some values
            handles.indexset = [];
            ClearPoints_Callback(hObject,[],handles);
            handles = guidata(hObject);
            handles.points = [];
            handles.normal = 1;
            handles.indexset = [];
            line = ['Drew edited region of object' num2str(boundaryindex)];
            PrintText(hObject,handles,line);
            guidata(hObject,handles);
            
            % In case of a second lesion, we put another listener in now
            set(handles.editedregion,'ButtonDownFcn',...
                {@SecondSelect,boundaryindex,handles});
        
        else
            % continue doing what we normally do when image is clicked
            ImageClickCallback(hObject,eventdata);
        end
        
    end
    
elseif (strcmp(handles.pointmode, 'brain'))
    coordinates = get(gca,'CurrentPoint');
    handles.brain = coordinates(1,1:2);
    plot(handles.brain(1,1),handles.brain(1,2),'x','LineWidth',3);
    handles.pointmode = 'midline';
    line = 'Exit Brain Select Mode - Enter Midline Draw Mode ';
    PrintText(hObject,handles,line);
    set(handles.LineModeText,'String','Midline');
    
elseif (strcmp(handles.pointmode, 'editing'))
    coordinates = get(gca,'CurrentPoint');
    clickedpoint =  coordinates(1,1:2);
    real = handles.real{object,1};
    dist = Inf;
    [m,~] = size(real);
    pointindex = 0;
    for i = 1:m
        check = mydistance(real(i,:),clickedpoint);
        if check < dist
            dist = check;
            pointindex = i;
        end
    end
    
    if (handles.edit == 0)
       handles.edit = pointindex;
       handles.points = cat(1,handles.points,clickedpoint);
       guidata(hObject,handles);
       line = 'First edit point selected ';
       PrintText(hObject,handles,line);
    elseif (pointindex > handles.edit)
        handles.real{object,1} = cat(1,real(1:handles.edit,:), ...
            handles.points,real(pointindex:end,:));
        
        r = handles.real{object,1};
        toss = handles.pHandle(object);
        handles.pHandle(object) = plot(r(:,1),r(:,2),'r','LineWidth',3);
        set(handles.pHandle(object),'ButtonDownFcn',{@PointSelect,object});
        handles.edit = 0;
        [m,~] = size(handles.points);
        if m > 1
            for i = 1:m-1
                delete(handles.segment{i})
            end
        end
        handles.points = [];
        handles.normal = 1;
        guidata(hObject,handles);
        delete(toss);
        line = 'New boundary drawn ';
        PrintText(hObject,handles,line);
    else
        handles.real{object,1} = cat(1,real(pointindex:handles.edit,:), ...
            handles.points,real(pointindex,:));
        r = handles.real{object,1};
        toss = handles.pHandle(object);
        handles.pHandle(object) = plot(r(:,1),r(:,2),'r','LineWidth',3);
        set(handles.pHandle(object),'ButtonDownFcn',{@PointSelect,object});
        handles.edit = 0;
        [m,~] = size(handles.points);
        if m > 1
            for i = 1:m-1
                delete(handles.segment{i});
            end
        end
        handles.points = [];
        handles.normal = 1;
        guidata(hObject,handles);
        delete(toss);
        line = 'New boundary drawn ';
        PrintText(hObject,handles,line);
    end
    guidata(hObject,handles);
else
    ImageClickCallback(hObject,eventdata); 
end

function SecondSelect (hObject,~, object, ~)
handles = guidata(hObject);

if ~(handles.mode) % in Area mode
    coordinates = get(gca,'CurrentPoint');
    clickedpoint = coordinates(1,1:2);

    % This time we want the boundset we made at the last edit
    boundset2 = handles.boundset{object,1};
    pointindex = indexOfPlot(boundset2,clickedpoint);

    % We only want 2 points on the boundary. More is taken as a new set
    [~, n] = size(handles.indexset);
    if n >= 2 
        handles.indexset = [];
    end

    % Now we add the point, provide feedback about where we are  
    handles.indexset = [handles.indexset,pointindex];
    handles.object = object; % To keep track of used bound
    set(handles.ObjectArea,'String',num2str(object));

    % Print feedback
    [~, n] = size(handles.indexset);
    line = ['Second Boundary ' num2str(n) ' selected on object ' ...
        num2str(object)];
    PrintText(hObject,handles,line);

    if n == 2 
        % Similar to PointSelect except with green line
        boundaryindex = handles.object;
        boundary = handles.boundset{boundaryindex,1};
        boundset = cat(1,boundary(1:handles.indexset(1,1),1:2),...
            handles.points, boundary(handles.indexset(1,2):end,1:2));
        handles.editedregion = plot(boundset(:,1),boundset(:,2),...
            'g','LineWidth',3);
        handles.boundset{boundaryindex,1} = boundset;

        % We re-initialize some values
        handles.indexset = [];
        handles.backup = handles.points;
        handles.points = [];
        handles.normal = 1;
        line = ['Drew edited region of object ' num2str(boundaryindex)];
        PrintText(hObject,handles,line);
        guidata(hObject,handles);
    else
        % continue doing what we normally do when image is clicked
        ImageClickCallback(hObject,[]);
    end
end

% --- Get click points on plot
function ImageClickCallback (hObject, ~)

handles = guidata(hObject);
axesHandle  = get(hObject,'Parent');
coordinates = get(axesHandle,'CurrentPoint'); 
coordinates = coordinates(1,1:2);
point = [coordinates(1) coordinates(2)]; 

%If we're in midline mode, we prepare to draw the midline
if (strcmp(handles.pointmode, 'midline'))
    handles.midlinepoints = cat(1,handles.midlinepoints,point);
    [m,~] = size(handles.midlinepoints);
    if (m == 2)
        p1 = handles.midlinepoints(1,:);
        p2 = handles.midlinepoints(2,:);
        yDiff = p2(2) - p1(2);
        xDiff = p2(1) - p1(1);
        m = yDiff/xDiff;
        x = -100000:1:100000;
        y = m*x - m*p1(1) + p1(2); % equation for line using 2 points
        handles.midline = plot(x,y,'k','LineWidth',3);
        handles.pointmode = 'normal';
        line = 'Exit Midline Draw Mode - Enter Normal Draw Mode ';
        PrintText(hObject,handles,line);
        set(handles.LineModeText,'String','Normal');
    end
    
%Otherwise, we draw regular points, which we keep track of to undo
elseif (strcmp(handles.pointmode, 'normal'))
    handles.points = cat(1,handles.points,point);
    [m,~] = size(handles.points);
    if m > 1
        i = handles.normal;
        x = [handles.points(i,1) handles.points(i+1,1)];
        y = [handles.points(i,2) handles.points(i+1,2)];
        handles.segment{i} = plot(x,y,'b','LineWidth',3);
        handles.normal = i + 1;
        
        % if we're in length mode, we only need 2 points to move on
        if(handles.mode)
            handles.pointmode = 'certain';
            line = 'Exit Normal Draw Mode - Enter Certain Draw Mode ';
            PrintText(hObject,handles,line);
            set(handles.LineModeText,'String','Certain');
        end
    end
    
%     % Show length of line
%     length = mydistance(handles.points(1,:),handles.points(end,:))...
%         *handles.cal;
%     line = ['Length of lesion is ' num2str(length * 1000) 'mm under ' ...
%         num2str(handles.cal*1000000) ' um/pixel'];
%     PrintText(hObject,handles,line);
    
% we do similar things to the certain lesion's points
elseif (strcmp(handles.pointmode, 'certain'))
    handles.cpoints = cat(1,handles.cpoints,point);
    [m,~] = size(handles.cpoints);
    if m > 1
        i = handles.certain;
        x = [handles.cpoints(i,1) handles.cpoints(i+1,1)];
        y = [handles.cpoints(i,2) handles.cpoints(i+1,2)];
        handles.csegment{i} = plot(x,y,'g','LineWidth',3);
        handles.certain = i + 1;
        
        % if we're in length mode, we only need 2 points to move on
        if(handles.mode)
            handles.pointmode = 'brain';
            line = 'Exit Certain Draw Mode - Normal Brain Select Entered ';
            PrintText(hObject,handles,line);
            set(handles.LineModeText,'String','Brain');
        end
    end
elseif (strcmp(handles.pointmode, 'editing'))
    handles.points = cat(1,handles.points,point);
    [m,~] = size(handles.points);
    if m > 1
        i = handles.normal;
        x = [handles.points(i,1) handles.points(i+1,1)];
        y = [handles.points(i,2) handles.points(i+1,2)];
        handles.segment{i} = plot(x,y,'b','LineWidth',3);
        handles.normal = i + 1; 
    end
end 
guidata(hObject,handles);

% Used to print text to the edit box under the plot. When the line is sent,
% the handles are needed to find the box and the line must be in the form
% of line = ['string' ] 
function PrintText(hObject,handles,line)
string = get(handles.DisplayInfo,'String');
result = [{line};string];
set(handles.DisplayInfo,'String',result);
guidata(hObject,handles);

% --- Executes on button press in AreaCall.
function AreaCall_Callback(hObject, ~, handles) %#ok<DEFNU>
% hObject    handle to AreaCall (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Find the object that will be plotted and highlight for checking
object = handles.object;
boundsetplot = handles.boundset{object,1};
realplot = handles.real{object,1};
plot(boundsetplot(:,1),boundsetplot(:,2),'g');

% Find the areas of the shapes and the difference is the lesion area 
editedArea = polyarea(boundsetplot(:,1),boundsetplot(:,2));
boundedArea = polyarea(realplot(:,1),realplot(:,2));
lesionArea = abs(editedArea - boundedArea);

%use the calibration to know the actual size
calArea = lesionArea * handles.cal^2;
calAreamm = calArea * 1000 * 1000;% convert from m^2 to mm^2
line =  ['Lesioned area of object ' num2str(object) ' equals: ' ...
    num2str(calAreamm) ' mm^2'];
PrintText(hObject,handles,line);
guidata(hObject,handles);


function ObjectArea_Callback(hObject, ~, handles) %#ok<DEFNU>
% hObject    handle to ObjectArea (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ObjectArea as text
%        str2double(get(hObject,'String')) returns contents of ObjectArea as a double
handles.object = str2double(get(hObject,'String'));
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function ObjectArea_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
% hObject    handle to ObjectArea (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in Import.
function Import_Callback(hObject, ~, handles)
% hObject    handle to Import (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Intialize/clear many values
handles.points = [];
handles.cpoints = [];
handles.midlinepoints = [];
handles.indexset = [];
handles.boundset = [];
handles.edit = 0;
handles.normal = 1;
handles.certain = 1;

if (handles.mode)
    handles.pointmode = 'brain';
    set(handles.LineModeText,'String','Brain');
end

cla reset %clear axes
filenameindex = handles.filenameindex; %find out what file we need
filename = strcat(handles.pathname,handles.filenames{filenameindex});

guidata(hObject,handles);
CalSet_Callback(hObject, [], handles)
handles = guidata(hObject);

% we can downsample the image right now to make it load faster. The
% variable idownsample means initial downsample
handles.image = imread(filename);
handles.imageHandle = imshow(handles.image);

% draw scale bar
hold on
scale =  0.001 / handles.cal; %scale in pixel/mm
y = ones([1,scale]) .* 50;
x = 50:1:scale+49;
plot(x,y,'c');
hold off

% now we set listeners on the image so we can respond to clicks
set(handles.imageHandle,'ButtonDownFcn', @ImageClickCallback);
line = ['Image of file ' handles.filenames{filenameindex} ' imported'];
PrintText(hObject,handles,line);
guidata(hObject,handles);


function Threshold_Callback(hObject, ~, handles) %#ok<DEFNU>
% hObject    handle to Threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Threshold as text
%        str2double(get(hObject,'String')) returns contents of Threshold as a double
handles.threshold = str2double(get(hObject,'String'));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function Threshold_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
% hObject    handle to Threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function DisplayInfo_Callback(~, ~, ~) %#ok<DEFNU>
% hObject    handle to DisplayInfo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DisplayInfo as text
%        str2double(get(hObject,'String')) returns contents of DisplayInfo as a double


% --- Executes during object creation, after setting all properties.
function DisplayInfo_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
% hObject    handle to DisplayInfo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in SelectFileButton.
function SelectFileButton_Callback(hObject, ~, handles) %#ok<DEFNU>
% hObject    handle to SelectFileButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% allows dialogue window to find mutliple files
[handles.filenames,handles.pathname,~] = uigetfile('*.tif','MultiSelect','on');

% if there is only one file, put it in a cell to conform to code
if ~iscellstr(handles.filenames)
    handles.filenames = {handles.filenames};
end

handles.filenameindex = 1;

% show the user what file we have
name = strtok(handles.filenames{1}, '_'); % extracts only slide name
set(handles.FileIndexTextBox,'String',name);
guidata(hObject,handles);

function FileIndexTextBox_Callback(hObject, ~, handles) %#ok<DEFNU>
% hObject    handle to FileIndexTextBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of FileIndexTextBox as text
%        str2double(get(hObject,'String')) returns contents of FileIndexTextBox as a double
handles.filenameindex = str2double(get(hObject,'String'));
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function FileIndexTextBox_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
% hObject    handle to FileIndexTextBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in NextIndexButton.
function NextIndexButton_Callback(hObject, ~, handles)
% hObject    handle to NextIndexButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[~,n] = size(handles.filenames);
if (n > handles.filenameindex)
    handles.filenameindex = handles.filenameindex + 1;
    filename = handles.filenames{handles.filenameindex};
    name = strtok(filename, '_'); % extracts only slide name
    set(handles.FileIndexTextBox,'String',name);
    guidata(hObject,handles);
    RowNext_Callback(hObject,[],handles);
else
    line = 'No further files in queue';
    PrintText(hObject,handles,line);
end
    
% --- Executes on button press in PreviousIndexButton.
function PreviousIndexButton_Callback(hObject, ~, handles) %#ok<DEFNU>
% hObject    handle to PreviousIndexButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if (handles.filenameindex > 1)
    handles.filenameindex = handles.filenameindex - 1;
    filename = handles.filenames{handles.filenameindex};
    name = strtok(filename, '_'); % extracts only slide name
    set(handles.FileIndexTextBox,'String',name);
    guidata(hObject,handles);
    RowBack_Callback(hObject,[],handles);
else
    line = 'No files before this point';
    PrintText(hObject,handles,line);
end

% --- Executes on button press in DownSample.
function DownSample_Callback(hObject, ~, handles) %#ok<DEFNU>
% hObject    handle to DownSample (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.image = imresize(handles.image,handles.downsample);
cla;
handles.cal = handles.cal / handles.downsample;
set(handles.CalText,'String',num2str(handles.cal*1000000));
imageHandle = imshow(handles.image);
set(imageHandle,'ButtonDownFcn', @ImageClickCallback);
line = ['Downsampled by ' num2str(handles.downsample)];
PrintText(hObject,handles,line);

% draw scale bar
hold on
scale =  0.001 / handles.cal; %scale in pixel/mm
y = ones([1,scale]) .* 50;
x = 50:1:scale+49;
plot(x,y,'c');
hold off

guidata(hObject,handles);

function DownSampleText_Callback(hObject, ~, handles) %#ok<DEFNU>
% hObject    handle to DownSampleText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DownSampleText as text
%        str2double(get(hObject,'String')) returns contents of DownSampleText as a double
handles.downsample = str2double(get(hObject,'String'));
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function DownSampleText_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
% hObject    handle to DownSampleText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function NoiseCancel_Callback(hObject, ~, handles) %#ok<DEFNU>
% hObject    handle to NoiseCancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of NoiseCancel as text
%        str2double(get(hObject,'String')) returns contents of NoiseCancel as a double
handles.noisecancel = str2double(get(hObject,'String'));
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function NoiseCancel_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
% hObject    handle to NoiseCancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function CalText_Callback(hObject, ~, handles) %#ok<DEFNU>
% hObject    handle to CalText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CalText as text
%        str2double(get(hObject,'String')) returns contents of CalText as a double
cal = str2double(get(hObject,'String')); %cal in um/pixel
handles.cal = cal / 1000000; %cal in m/pixel
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function CalText_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
% hObject    handle to CalText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ToggleLineMode.
function ToggleLineMode_Callback(hObject, ~, handles) %#ok<DEFNU>
% hObject    handle to ToggleLineMode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Order of toggle: normal > certain > midline > normal .....
if (strcmp(handles.pointmode, 'certain'))
    handles.pointmode = 'brain';
    line = 'Exit Certain Draw Mode - Normal Brain Select Entered ';
    PrintText(hObject,handles,line);
    set(handles.LineModeText,'String','Brain');
    
elseif (strcmp(handles.pointmode, 'brain'))
    handles.pointmode = 'midline';
    line = 'Exit Brain Select Mode - Enter Midline Draw Mode ';
    PrintText(hObject,handles,line);
    set(handles.LineModeText,'String','Midline');
    
elseif (strcmp(handles.pointmode, 'midline'))
    handles.pointmode = 'normal';
    line = 'Exit Midline Draw Mode - Enter Normal Draw Mode ';
    PrintText(hObject,handles,line);
    set(handles.LineModeText,'String','Normal');

else
    handles.pointmode = 'certain';
    line = 'Exit Normal Draw Mode - Enter Certain Draw Mode ';
    PrintText(hObject,handles,line);
    set(handles.LineModeText,'String','Certain');
end
guidata(hObject,handles);

% --- Executes on button press in UndoPoint.
function UndoPoint_Callback(hObject, ~, handles) %#ok<DEFNU>
% hObject    handle to UndoPoint (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Actual purpose below
i = handles.normal - 1;
if i > 0
    delete(handles.segment{i});
    handles.points = handles.points(1:i,1:2);
    handles.normal = i;
else
    handles.points = [];
    handles.normal = 1;
    handles.indexset = [];
    line = 'Boundary cleared';
    PrintText(hObject,handles,line);
end
guidata(hObject,handles);


% --- Executes on button press in UndoEditedDraw.
function UndoEditedDraw_Callback(hObject, ~, handles) %#ok<DEFNU>
% hObject    handle to UndoEditedDraw (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
delete(handles.editedregion);
guidata(hObject,handles);



% --- Executes on button press in ClearPoints.
function ClearPoints_Callback(hObject, ~, handles)
% hObject    handle to ClearPoints (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if (strcmp(handles.pointmode, 'normal'))
    i = handles.normal - 1;
    for l = 1:i
        delete(handles.segment{l});
    end
    handles.points = [];
    handles.normal = 1;

elseif (strcmp(handles.pointmode, 'midline'))
    delete(handles.midline);
    handles.midlinepoints = [];

elseif (strcmp(handles.pointmode, 'certain'))
    i = handles.certain - 1;
    for l = 1:i
        delete(handles.csegment{l});
    end
    handles.cpoints = [];
    handles.certain = 1;
end
guidata(hObject,handles);

% --- Executes on button press in CropButton.
function CropButton_Callback(hObject, ~, handles) %#ok<DEFNU>
% hObject    handle to CropButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.image = imcrop;
delete(handles.imageHandle);
handles.imageHandle = imshow(handles.image);

% draw scale bar
hold on
scale =  0.001 / handles.cal; %scale in pixel/mm
y = ones([1,scale]) .* 50;
x = 50:1:scale+49;
plot(x,y,'c');
hold off

% now we set listeners on the image so we can respond to clicks
set(handles.imageHandle,'ButtonDownFcn', @ImageClickCallback);
line = ['Image of file ' handles.filenames{handles.filenameindex} ' cropped'];
PrintText(hObject,handles,line);
guidata(hObject,handles);

% --- Executes on button press in ShowBW.
function ShowBW_Callback(~, ~, handles) %#ok<DEFNU>
% hObject    handle to ShowBW (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
I = handles.image;
BW = im2bw(I,handles.threshold);
BW = ~BW; %Invert because samples are labeled as background accidentally
figure, imshow(BW);

function AppendPlot(hObject,handles)
% Keep track of what slice you're on with row. After storing the
% information in the proper location we move onto the next slide and
% therefore the next file we imported.
row = handles.row;

if (handles.mode) % actions in length mode

    % Find the distance between the lesion and the midline and length of lesion
    [m,~] = size(handles.points);
    if ~(m == 0)
        % find distance from points on line to the midpoint
        d = distancetoline(handles.points(end,:),handles.midlinepoints(1,:),handles.midlinepoints(2,:));
        d2 = distancetoline(handles.points(1,:),handles.midlinepoints(1,:),handles.midlinepoints(2,:));

        % shorter distance is the distance to lesion
        if (d > d2)
            dist = d2;
            len = d;
        else 
            dist = d;
            len = d2;
        end

        % adjust to from pixel to meters
        distance = dist * handles.cal;
        line = ['Lesion is ' num2str(distance*1000) 'mm from the midline'];
        PrintText(hObject,handles,line);

        % Find length of lesion area and convert to meters
        length = mydistance(handles.points(1,:),handles.points(end,:))*handles.cal;
        line = ['Length of lesion is ' num2str(length * 1000) 'mm under ' num2str(handles.cal*1000000) ' um/pixel'];
        PrintText(hObject,handles,line);
        
        % Find length of lesion area and convert to meters
        length = (len - dist)*handles.cal;
        line = ['Length of lesion is ' num2str(length * 1000) 'mm under ' num2str(handles.cal*1000000) ' um/pixel'];
        PrintText(hObject,handles,line);

        cdistance = distance;
        clength = length;

        [cm,~] = size(handles.cpoints);
        if ~(cm == 0)
            % find how far the needed end of the certain lesion is from the midline
            cd = distancetoline(handles.cpoints(end,:),handles.midlinepoints(1,:),handles.midlinepoints(2,:));
            cd2 = distancetoline(handles.cpoints(1,:),handles.midlinepoints(1,:),handles.midlinepoints(2,:));

            % shorter distance is the distance to lesion
            if (cd > cd2)
                dist = cd2;
                len = cd;
            else 
                dist = cd;
                len = cd2;
            end
            
            cdistance = dist * handles.cal;
            line = ['The certain lesion is ' num2str(cdistance*1000) 'mm from the midline'];
            PrintText(hObject,handles,line);

            % Find length of the certain lesion area
            clength = (len - dist)*handles.cal;
            line = ['Length of certain lesion is ' num2str(clength * 1000) 'mm under ' num2str(handles.cal*1000000) ' um/pixel'];
            PrintText(hObject,handles,line);
        end

    % if nothing was inputed, then we have 0s for values
    else 
        distance = 0;
        cdistance = 0;
        length = 0;
        clength = 0;
    end

    % find how far the end of the brain is from the midline
    b = distancetoline(handles.brain,handles.midlinepoints(1,:),handles.midlinepoints(2,:));
    brain = b * handles.cal;
    line = ['The brain ends ' num2str(brain*1000) 'mm from the midline'];
    PrintText(hObject,handles,line);

    % We store the information with 10000x to give them units of 100um so we
    % have a finer degree of information since indexing use whole numbers
    handles.colorPlots(row,1) = distance * 10000;
    handles.colorPlots(row,2) = length * 10000;
    handles.colorPlots(row,3) = cdistance * 10000;
    handles.colorPlots(row,4) = clength * 10000;
    handles.colorPlots(row,5) = brain * 10000;
    handles.colorPlots(row,6) = 0;
    % Mark Bregma's spot
    if (row == handles.bregma)
        handles.colorPlots(row,6) = 1;
    end
    
    line = ['Appended information to row ' num2str(row)];
    PrintText(hObject,handles,line);
    guidata(hObject,handles);
    
else % otherwise we're in area mode so we have a separate set of actions

    % We store the information with 10000x to give them units of 100um so we
    % have a finer degree of information since indexing use whole numbers
    handles.plot{row,1} = handles.real;
    handles.plot{row,2} = handles.boundset;
    handles.plot{row,3} = handles.cal;

    line = ['Appended information to row ' num2str(row)];
    PrintText(hObject,handles,line);
    guidata(hObject,handles);
end

% --- Executes on button press in Append.
function Append_Callback(hObject, ~, handles) %#ok<DEFNU>
% hObject    handle to Append (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

AppendPlot(hObject,handles);
handles = guidata(hObject);
handles.row = mod(handles.row,100);
% Move to next image
NextIndexButton_Callback(hObject, [], handles)
handles = guidata(hObject);
Import_Callback(hObject, [], handles)

% --- Executes on button press in Save.
function Save_Callback(~, ~, handles) %#ok<DEFNU>
% hObject    handle to Save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if (handles.mode)
    colorplot = handles.colorPlots;
    save(strcat('ColorPlots\',handles.outFile),'colorplot');
    MakeColorPlot(colorplot,handles.outFile);
else 
    plots = handles.plot;
    save(strcat('CartoonPlots\',handles.outFile),'plots');
    ShowLesions(plots,handles.outFile);
end

function Filename_Callback(hObject, ~, handles) %#ok<DEFNU>
% hObject    handle to Filename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Filename as text
%        str2double(get(hObject,'String')) returns contents of Filename as a double
handles.outFile = get(hObject,'String');
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function Filename_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
% hObject    handle to Filename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function LineModeText_Callback(~, ~, ~) %#ok<DEFNU>
% hObject    handle to LineModeText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of LineModeText as text
%        str2double(get(hObject,'String')) returns contents of LineModeText as a double


% --- Executes during object creation, after setting all properties.
function LineModeText_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
% hObject    handle to LineModeText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ModeText_Callback(~, ~, ~) %#ok<DEFNU>
% hObject    handle to ModeText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ModeText as text
%        str2double(get(hObject,'String')) returns contents of ModeText as a double

% --- Executes during object creation, after setting all properties.
function ModeText_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
% hObject    handle to ModeText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Toggle. Switch between Length and Area
function Toggle_Callback(hObject, ~, handles)
% hObject    handle to Toggle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if (handles.mode)
    handles.mode = 0;
    set(handles.ModeText,'String','Area');
    handles.pointmode = 'normal';
    set(handles.LineModeText,'String','Normal');
else 
    handles.mode = 1;
    set(handles.ModeText,'String','Length');
    handles.pointmode = 'brain';
    set(handles.LineModeText,'String','Brain');
end
guidata(hObject,handles);


% --- Executes on button press in EditBoundary.
function EditBoundary_Callback(hObject, ~, handles) %#ok<DEFNU>
% hObject    handle to EditBoundary (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if (strcmp(handles.pointmode,'editing'))
    if handles.mode
        handles.pointmode = 'brain';
        line = 'Enter brain mode ';
        PrintText(hObject,handles,line);
        set(handles.LineModeText,'String','Brain');
        guidata(hObject,handles);
    else
        handles.pointmode = 'normal';
        line = 'Enter Lesion mode ';
        PrintText(hObject,handles,line);
        set(handles.LineModeText,'String','Normal');
        guidata(hObject,handles);
    end
else
    handles.pointmode = 'editing';
    line = 'Edit boundaries now ';
    PrintText(hObject,handles,line);
    set(handles.LineModeText,'String','Editing');
    guidata(hObject,handles);
end


% --- Executes on button press in UndoBoundary.
function UndoBoundary_Callback(hObject, ~, handles) %#ok<DEFNU>
% hObject    handle to UndoBoundary (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[~,n] = size(handles.pHandle);
for i = 1:n
    delete(handles.pHandle(i))
end
guidata(hObject,handles);

% --- Executes on button press in ImportOld.
function ImportOld_Callback(hObject, ~, handles) %#ok<DEFNU>
% hObject    handle to ImportOld (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Select file and load
[filename,pathname,~] = uigetfile('*.mat');
location = strcat(pathname,filename);
file = load(location);
handles.outFile = strrep(filename,'.mat','');
set(handles.Filename,'String',handles.outFile);

% Select the correct mode
if isfield(file,'colorplot') && ~(handles.mode) || isfield(file,'plots') && (handles.mode)
    Toggle_Callback(hObject, [], handles);
end
handles = guidata(hObject);

% Clear previous data
handles.colorPlots = [];
handles.plot = [];

% If we're in length mode, we have processed and results files
if (handles.mode)
    handles.colorPlots = file.colorplot;   
    line = ['Loaded ' filename ' in length mode'];
    PrintText(hObject,handles,line);
    figure,MakeColorPlot(handles.colorPlots,handles.outFile);
    
% Otherwise we are in area mode which contains sets of plots
else
    handles.plot = file.plots;
    line = ['Loaded ' filename ' in area mode'];
    PrintText(hObject,handles,line);
    ShowLesions(handles.plot,handles.outFile);
end

guidata(hObject,handles);

% --- Executes on button press in RowNext.
function RowNext_Callback(hObject, ~, handles)
% hObject    handle to RowNext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.row = handles.row + 1;
set(handles.Row,'String',num2str(handles.row));
guidata(hObject,handles);

% --- Executes on button press in RowBack.
function RowBack_Callback(hObject, ~, handles)
% hObject    handle to RowBack (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.row = handles.row - 1;
set(handles.Row,'String',num2str(handles.row));
guidata(hObject,handles);


function Row_Callback(hObject, ~, handles) %#ok<DEFNU>
% hObject    handle to Row (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Row as text
%        str2double(get(hObject,'String')) returns contents of Row as a double
handles.row = str2double(get(hObject,'String'));
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function Row_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
% hObject    handle to Row (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function AutoCal_Callback(~, ~, ~) %#ok<DEFNU>
% hObject    handle to AutoCal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of AutoCal as text
%        str2double(get(hObject,'String')) returns contents of AutoCal as a double


% --- Executes during object creation, after setting all properties.
function AutoCal_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
% hObject    handle to AutoCal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in CalSet.
function CalSet_Callback(hObject, ~, handles)
% hObject    handle to CalSet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
filenameindex = handles.filenameindex; %find out what file we need
filename = handles.filenames{filenameindex};

[~,remain] = strtok(filename, '_'); % removes slide name

% Auto set calibration
if ~isempty(strfind(filename,'cal'))
    [~,remain] = strtok(remain, '_'); % removes 'cal'
    [cal,remain] = strtok(remain,'_'); % gets start of calibration
    cal = strrep(cal,'.tif','');
    handles.cal = str2double(cal) / 1000000;
    set(handles.CalText,'String',num2str(handles.cal*1000000));
    set(handles.AutoCal,'String','1');
end

% If the file has been changed by LabScale, we must adjust calibration
if ~isempty(strfind(filename,'down'))
    [~,remain] = strtok(remain,'_'); % removes 'down'
    autocal = strtok(remain,'_'); % saves downsample rate
    autocal = strrep(autocal,'.tif','');
    set(handles.AutoCal,'String',autocal);
    handles.cal = handles.cal / str2double(autocal);
    set(handles.CalText,'String',num2str(handles.cal*1000000));
end

guidata(hObject,handles);


% --- Executes on button press in Second.
function Second_Callback(hObject, ~, handles) %#ok<DEFNU>
% hObject    handle to Second (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% To allow mutliple lesions on one brain, we append the information from
% the first and start collecting information for the second in a new row 
AppendPlot(hObject,handles);
handles = guidata(hObject);
handles.row = handles.row + 100;

% reset variables depending on the mode
if (handles.mode) % in Length mode
    handles.brain = 0;
    delete(handles.midline); % remove to allow for new one to be drawn
    handles.midlinepoints = [];
    handles.points = [];
    handles.cpoints = [];
    handles.normal = 1;
    handles.certain = 1;  
end

guidata(hObject,handles);


% --- Executes on button press in Bregma.
function Bregma_Callback(hObject, ~, handles) %#ok<DEFNU>
% hObject    handle to Bregma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.bregma = handles.row;
guidata(hObject,handles);