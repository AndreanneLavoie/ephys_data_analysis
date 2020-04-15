function varargout = electrodeSelectorSingle(varargin)
% ELECTRODESELECTORSINGLE MATLAB code for electrodeSelectorSingle.fig
%      ELECTRODESELECTORSINGLE, by itself, creates a new ELECTRODESELECTORSINGLE or raises the existing
%      singleton*.
%
%      H = ELECTRODESELECTORSINGLE returns the handle to a new ELECTRODESELECTORSINGLE or the handle to
%      the existing singleton*.
%
%      ELECTRODESELECTORSINGLE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ELECTRODESELECTORSINGLE.M with the given input arguments.
%
%      ELECTRODESELECTORSINGLE('Property','Value',...) creates a new ELECTRODESELECTORSINGLE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before electrodeSelectorSingle_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to electrodeSelectorSingle_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help electrodeSelectorSingle

% Last Modified by GUIDE v2.5 29-Aug-2019 17:37:55

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @electrodeSelectorSingle_OpeningFcn, ...
                   'gui_OutputFcn',  @electrodeSelectorSingle_OutputFcn, ...
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


% --- Executes just before electrodeSelectorSingle is made visible.
function electrodeSelectorSingle_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to electrodeSelectorSingle (see VARARGIN)

% Choose default command line output for electrodeSelectorSingle
handles.output = hObject;
contents = cellstr(get(handles.popupmenu1,'String'));
handles.ampA = contents{get(handles.popupmenu1,'Value')};
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes electrodeSelectorSingle wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = electrodeSelectorSingle_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.ampA;
varargout{2} = handles.okstatus;

delete(handles.figure1);




% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

contents = cellstr(get(hObject,'String')); %returns popupmenu1 contents as cell array
ampA = contents{get(hObject,'Value')}; %returns selected item from popupmenu1
handles.ampA = ampA;
% Update handles structure
guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in okbutton.
function okbutton_Callback(hObject, eventdata, handles)
% hObject    handle to okbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.okstatus = 1;
% Update handles structure
guidata(hObject, handles);

% Resume
hFig = ancestor(hObject,'Figure');
if isequal(get(hFig, 'waitstatus'), 'waiting')
    % The GUI is still in UIWAIT
    uiresume(handles.figure1);
end


% --- Executes on button press in cancelbutton.
function cancelbutton_Callback(hObject, eventdata, handles)
% hObject    handle to cancelbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.okstatus = 0;
% Update handles structure
guidata(hObject, handles);

% Resume
hFig = ancestor(hObject,'Figure');
if isequal(get(hFig, 'waitstatus'), 'waiting')
    % The GUI is still in UIWAIT
    uiresume(handles.figure1);
end


% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function figure1_WindowButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object deletion, before destroying properties.
function Untitled_1_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on key release with focus on figure1 or any of its controls.
function figure1_WindowKeyReleaseFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.FIGURE)
%	Key: name of the key that was released, in lower case
%	Character: character interpretation of the key(s) that was released
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) released
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on key press with focus on figure1 and none of its controls.
function figure1_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
