function varargout = citaToolbox(varargin)
% CITATOOLBOX MATLAB code for citaToolbox.fig
%      CITATOOLBOX, by itself, creates a new CITATOOLBOX or raises the existing
%      singleton*.
%
%      H = CITATOOLBOX returns the handle to a new CITATOOLBOX or the handle to
%      the existing singleton*.
%
%      CITATOOLBOX('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CITATOOLBOX.M with the given input arguments.
%
%      CITATOOLBOX('Property','Value',...) creates a new CITATOOLBOX or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before citaToolbox_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to citaToolbox_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help citaToolbox

% Last Modified by GUIDE v2.5 11-Jul-2019 11:31:07

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @citaToolbox_OpeningFcn, ...
    'gui_OutputFcn',  @citaToolbox_OutputFcn, ...
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


% --- Executes just before citaToolbox is made visible.
function citaToolbox_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to citaToolbox (see VARARGIN)

%Display all handles fields (dev utility)
disp(handles);

%Hide axes (use axes as image containers)
set(handles.citaLogo,'visible','off');
set(handles.mooreLogo,'visible','off');
set(handles.riverImage,'visible','off');

%Disable buttons until image is selected
set(handles.reclassifyApplyButton,'Enable','off');
set(handles.reclassifySaveButton,'Enable','off');
set(handles.sizeFilterApplyButton,'Enable','off');
set(handles.sizeFilterSaveButton,'Enable','off');
set(handles.polygonizeButton,'Enable','off');
set(handles.polygonizeSaveButton,'Enable','off');
set(handles.centerlineButton,'Enable','off');
set(handles.centerlineSaveButton,'Enable','off');

%Set label and textfield initial values
set(handles.selectedImageText,'String','');
set(handles.reclassifyValue,'String','0.0');
set(handles.sizeFilterValue,'String','0.0');
set(handles.uploadShapeText,'String','');

%Read CITA and Moore logos images
[logoCita, ~, citaAlphaChannel] = imread('logoCita.png');
[logoMoore, ~, mooreAlphaChannel] = imread('logoMoore.png');
[imagePlaceholder, ~, placeholderAlphaChannel] = imread('logoCitaLarge.png');

%Add images to axes (alpha channel is used to show transparent images correctly)
axes(handles.citaLogo);
set(imshow(logoCita), 'AlphaData', citaAlphaChannel);

axes(handles.mooreLogo);
set(imshow(logoMoore), 'AlphaData', mooreAlphaChannel);

axes(handles.riverImage);
set(imshow(imagePlaceholder), 'AlphaData', placeholderAlphaChannel);

%Choose default command line output for citaToolbox
handles.output = hObject;

%Update handles structure
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = citaToolbox_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in selectImageButton.
function selectImageButton_Callback(hObject, eventdata, handles)
% hObject    handle to selectImageButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%File picker shows to select image (only format available for now is tif)
[fileName, pathName] = uigetfile({'*.tif';});
handles.pathName = pathName;
handles.fileName = fileName;

%Concat full path (file path + file name)
fullPath = strcat(pathName,fileName);
handles.fullPath = fullPath;

%User press cancel
if isequal(fullPath,0)
    %No action
else
    %imageData: double type matrix with index values
    %refmat: Referencing matrix
    %bbox: Bounding box
    %R: spatial referencing object
    [imageData, refmat, bbox] = geotiffread(fullPath);
    [~, R] = geotiffread(fullPath);
    
    %Mask outliers values from the matrix (some QGIS tif exports have out
    %of range values)
    imageData(imageData == 0) = 1;
    imageData(imageData >= 1) = 1;
    imageData(imageData < -1) = 1;
    imageData(isnan(imageData)) = 1;
    
    %save variables to handles
    handles.imageData = imageData;
    handles.refmat = refmat;
    handles.bbox = bbox;
    handles.R = R;
    
    %Set tif image to axes
    axes(handles.riverImage);
    imshow(imageData);
    
    %Set image name in label
    set(handles.selectedImageText, 'String', fileName);
    
    %Enable reclassify buttons and disable all the others
    set(handles.reclassifyApplyButton,'Enable','on');
    set(handles.reclassifySaveButton,'Enable','on');
    set(handles.sizeFilterApplyButton,'Enable','off');
    set(handles.sizeFilterSaveButton,'Enable','off');
    set(handles.polygonizeButton,'Enable','off');
    set(handles.polygonizeSaveButton,'Enable','off');
    set(handles.centerlineButton,'Enable','off');
    set(handles.centerlineSaveButton,'Enable','off');
    
    %save handles
    guidata(hObject, handles)
    
end


function reclassifyValue_Callback(hObject, eventdata, handles)
% hObject    handle to reclassifyValue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Get reclassify value entered by the user and save it to handles
reclassifyValue = get(hObject, 'String');
handles.reclassifyValue = reclassifyValue;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function reclassifyValue_CreateFcn(hObject, eventdata, handles)
% hObject    handle to reclassifyValue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in reclassifyApplyButton.
function reclassifyApplyButton_Callback(hObject, eventdata, handles)
% hObject    handle to reclassifyApplyButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%If no value entered, error message shows
if isempty(handles.reclassifyValue)
    warndlg("Enter numeric value", "Reclassify");
else
    %Get reclassify input and convert to string
    reclassifyThresholdDouble = str2double(handles.reclassifyValue);
    
    %If value is not numeric, error message shows
    if(isnan(reclassifyThresholdDouble))
        warndlg("Enter numeric value", "Reclassify");
    %If valid value is entered reclassification applies: pixels with lower value
    %than threshold returns a 0 value and pixels with higher values return 1.
    else
        %Temp variable init (copy of imageData in order to apply reclassification)
        reclassifiedData = handles.imageData;
        
        %RECLASSIFICATION
        
        %NDVI
        reclassifiedData(reclassifiedData <= reclassifyThresholdDouble) = 0;
        reclassifiedData(reclassifiedData > reclassifyThresholdDouble) = 1;
        reclassifiedData = ~reclassifiedData;
        
        %Save temp variable in handles
        handles.reclassifiedData = reclassifiedData;
        guidata(hObject, handles);
        
        %Show reclassified data as image
        axes(handles.riverImage);
        imshow(reclassifiedData);
        
        %Enable size filter buttons
        set(handles.sizeFilterApplyButton,'Enable','on');
        set(handles.sizeFilterSaveButton,'Enable','on');

    end
end


% --- Executes on button press in reclassifySaveButton.
function reclassifySaveButton_Callback(hObject, eventdata, handles)
% hObject    handle to reclassifySaveButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Open save file window (only tif file format for now)
[file,path] = uiputfile('*.tif', 'Save reclassified image', 'reclassify.tif');

if(file == 0)
    %No action
else
    geotiffwrite(strcat(path, file), handles.reclassifiedData, handles.R);
end


function sizeFilterValue_Callback(hObject, eventdata, handles)
% hObject    handle to sizeFilterValue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Get size filter value entered by the user and save it to handles
sizeFilterValue = get(hObject, 'String');
handles.sizeFilterValue = sizeFilterValue;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function sizeFilterValue_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sizeFilterValue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in sizeFilterApplyButton.
function sizeFilterApplyButton_Callback(hObject, eventdata, handles)
% hObject    handle to sizeFilterApplyButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%If not value entered, error messages shows
if isempty(handles.sizeFilterValue)
    warndlg("Enter numeric value (positive integer)", "Size filter");
else
    sizeFilterThresholdDouble = str2double(handles.sizeFilterValue);
    %If not numeric value entered, alert message shows
    if(isnan(sizeFilterThresholdDouble) || sizeFilterThresholdDouble < 0)
        warndlg("Enter numeric value (positive integer)", "Size filter");
    %If valid value is entered, connected pixels with a size value (as a
    %group) less than the value entered are filtered.
    else
        sizeFilteredData = bwareaopen(handles.reclassifiedData, sizeFilterThresholdDouble);
        
        %save size filtered data variable in handles
        handles.sizeFilteredData = sizeFilteredData;
        guidata(hObject, handles);
        
        %Show filtered image in axes
        axes(handles.riverImage);
        imshow(sizeFilteredData);
        
        %Enable polygonize buttons
        set(handles.polygonizeButton,'Enable','on');
        set(handles.polygonizeSaveButton,'Enable','on');

    end
end


% --- Executes on button press in sizeFilterSaveButton.
function sizeFilterSaveButton_Callback(hObject, eventdata, handles)
% hObject    handle to sizeFilterSaveButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Open save file window (only tif file format for now)
[file,path] = uiputfile('*.tif', 'Save size filtered image', 'sizeFiltered.tif');

if(file == 0)
    %No action
else
    geotiffwrite(strcat(path, file), handles.sizeFilteredData, handles.R);
end


% --- Executes on button press in polygonizeButton.
function polygonizeButton_Callback(hObject, eventdata, handles)
% hObject    handle to polygonizeButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Get the outline of the image (previously reclassified and size filtered)
outline = bwmorph(handles.sizeFilteredData,'remove');

%Get the location of the 1s in the matrix (row & column)
[rowNonZeroOutline, colNonZeroOutline] = find(outline);
c1Outline = rowNonZeroOutline(1);
r1Outline = colNonZeroOutline(1);

%Rearrange the points clockwise and going East
contourOutline = bwtraceboundary(outline,[c1Outline r1Outline],'E');
contourOutlineSize = size(contourOutline);

%Use reference matrix (refmat) to add geographic coordinates to the vector
%Create vector in handles to store the rearranged points
georeferencedContourOutline = zeros(1,2);

for i = 1:contourOutlineSize(1)
    latLongContourOutline = horzcat([contourOutline(i,1) contourOutline(i,2)], 1) * handles.refmat;
    georeferencedContourOutline = vertcat(georeferencedContourOutline, latLongContourOutline);
end

georeferencedContourOutline = georeferencedContourOutline(2:end,:);

%Save variables to handles
handles.georeferencedContourOutline = georeferencedContourOutline;
guidata(hObject, handles);
axes(handles.riverImage);

%Enable centerline buttons
set(handles.centerlineButton,'Enable','on');
set(handles.centerlineSaveButton,'Enable','on');

% plot(flipud(handles.georeferencedContourOutline(:,1)), flipud(handles.georeferencedContourOutline(:,2)),'Color', [0, 0, 0],'LineWidth',2);
imshow(outline);


% --- Executes on button press in polygonizeSaveButton.
function polygonizeSaveButton_Callback(hObject, eventdata, handles)
% hObject    handle to polygonizeSaveButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Open save file window (only shp file format for now)
[file,path] = uiputfile('*.shp', 'Save polygon', 'polygon.shp');

if(file == 0)
    %No action
else
    %Smooth data matrix
    handles.georeferencedContourOutline = sgolayfilt(handles.georeferencedContourOutline, 3, 5);
    
    %Create geovector with coordinates as vector (output file of type LINE)
    geoPointVectorContourOutline = geoshape(handles.georeferencedContourOutline(:,2), handles.georeferencedContourOutline(:,1));
    geoPointVectorContourOutline.Geometry = 'line';
    shapewrite(geoPointVectorContourOutline, strcat(path, file));
end


% --- Executes on button press in centerlineButton.
function centerlineButton_Callback(hObject, eventdata, handles)
% hObject    handle to centerlineButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Get the centerline (skeleton) of the image (previously reclassified and size filtered)
centerline = bwmorph(handles.sizeFilteredData,'skel',Inf);
%centerlineOnes = length(find(centerline));

S = bwmorph(handles.sizeFilteredData,'skel',Inf);
A = bwmorph(S,'branchpoints');
B = bwmorph(S,'endpoints');
C = imfuse(A, B);

%Flag to use code to remove spurs (still testing)
removeSpurs = 0;

if(removeSpurs == 1)
    
    I = centerline;

    %Alternative splitting method to 'branchpoint'
    %Use convolution to identify points with more than 2 neighboring pixels
    filter = [1 1 1;
              1 0 1;
              1 1 1];

    I_disconnect = I & ~(I & conv2(double(I), filter, 'same') > 2);

    cc = bwconncomp(I_disconnect);
    numPixels = cellfun(@numel,cc.PixelIdxList);
    [sorted_px, ind] = sort(numPixels);

    %Remove components shorter than threshold
    threshold  = 50;
    for ii=ind(sorted_px<threshold)
        cur_comp = cc.PixelIdxList{ii};
        I(cur_comp) = 0; 

        %Before removing component, check whether image is still connected
        full_cc = bwconncomp(I);
        if full_cc.NumObjects>1
            I(cur_comp) = 1; 
        end
    end

    %Clean up left over spurs
    I = bwmorph(I, 'spur');

    centerline = I;
    
end

[onesCenterlineX, onesCenterlineY] = find(centerline');
sizeOnesCenterline = size(onesCenterlineX);

georeferencedContour = zeros(1,2);

for i = 1:sizeOnesCenterline(1)
    latLongContour = horzcat([onesCenterlineY(i) onesCenterlineX(i)], 1) * handles.refmat;
    georeferencedContour = vertcat(georeferencedContour, latLongContour);
end

%Get the location of the 1s in the matrix (row & column)
%[rowNonZero, colNonZero] = find(centerline);

%c1 = rowNonZero(end);
%r1 = colNonZero(end);

%Rearrange the points clockwise and going East
%contour = bwtraceboundary(centerline, [c1 r1], 'E', 8, Inf, 'clockwise');
%contourSize = size(contour);
%assignin('base','contour',contour);


%Use reference matrix (refmat) to add geographic coordinates to the vector
%Create vector in handles to store the rearranged points
% georeferencedContour = zeros(1,2);
% 
% for i = 1:contourSize(1)
%     latLongContour = horzcat([contour(i,1) contour(i,2)], 1) * handles.refmat;
%     georeferencedContour = vertcat(georeferencedContour, latLongContour);
% end

georeferencedContour = georeferencedContour(2:end,:);
assignin('base','georeferencedContour',georeferencedContour);

%georeferencedContour = unique(georeferencedContour, 'stable', 'rows');

%Get repeated items index
% [~, ind] = unique(georeferencedContour(:, :), 'stable', 'rows');
% duplicate_ind = setdiff(1:size(georeferencedContour, 1), ind);

%Save variables to handles
handles.georeferencedContour = georeferencedContour;
%handles.georeferencedContour = sgolayfilt(handles.georeferencedContour, 3, 5);
guidata(hObject, handles);

%plot(georeferencedContour(:,1), georeferencedContour(:,2), 'Color', [0, 0, 0], 'LineWidth',2);
%scatter(georeferencedContour(:,1), georeferencedContour(:,2));
imshow(centerline);


% --- Executes on button press in centerlineSaveButton.
function centerlineSaveButton_Callback(hObject, eventdata, handles)
% hObject    handle to centerlineSaveButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Open save file window (only shp file format for now)
[file,path] = uiputfile('*.shp', 'Save centerline', 'centerline.shp');

if(file == 0)
    %No action
else
    %Create geovector with coordinates as vector (output file of type Point)
    geoPointVectorCenterline = geoshape(handles.georeferencedContour(:,2), handles.georeferencedContour(:,1));
    geoPointVectorCenterline.Geometry = 'point';
    %plot(handles.georeferencedContour(:,2), handles.georeferencedContour(:,1));
    shapewrite(geoPointVectorCenterline, strcat(path, file));
end


% --- Executes on button press in uploadButton.
function uploadButton_Callback(hObject, eventdata, handles)
% hObject    handle to uploadButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[file, path] = uigetfile({'*.shp';});

%User press cancel
if (file == 0)
    %No action
else
    S = shaperead(strcat(path, file));
    sizeS = size(S);
    for i = 1:sizeS(1)
        plot(S(i).X, S(i).Y, 'Color', [0, 0, 0]);
        hold on
    end
    hold off
    set(handles.uploadShapeText, 'String', file);
    
%   assignin('base','S',S);
%   display(S);
    %plot(S.X, S.Y);
  
end


% --- Executes on button press in test_button.
function test_button_Callback(hObject, eventdata, handles)
% hObject    handle to test_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
