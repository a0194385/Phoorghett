function varargout = UltrasoundImagingSimulator_Zheheng(varargin)
% ULTRASOUNDIMAGINGSIMULATOR_ZHEHENG MATLAB code for UltrasoundImagingSimulator_Zheheng.fig
%      ULTRASOUNDIMAGINGSIMULATOR_ZHEHENG, by itself, creates a new ULTRASOUNDIMAGINGSIMULATOR_ZHEHENG or raises the existing
%      singleton*.
%
%      H = ULTRASOUNDIMAGINGSIMULATOR_ZHEHENG returns the handle to a new ULTRASOUNDIMAGINGSIMULATOR_ZHEHENG or the handle to
%      the existing singleton*.
%
%      ULTRASOUNDIMAGINGSIMULATOR_ZHEHENG('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ULTRASOUNDIMAGINGSIMULATOR_ZHEHENG.M with the given input arguments.
%
%      ULTRASOUNDIMAGINGSIMULATOR_ZHEHENG('Property','Value',...) creates a new ULTRASOUNDIMAGINGSIMULATOR_ZHEHENG or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before UltrasoundImagingSimulator_Zheheng_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to UltrasoundImagingSimulator_Zheheng_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help UltrasoundImagingSimulator_Zheheng

% Last Modified by GUIDE v2.5 16-Sep-2020 11:39:14

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @UltrasoundImagingSimulator_Zheheng_OpeningFcn, ...
                   'gui_OutputFcn',  @UltrasoundImagingSimulator_Zheheng_OutputFcn, ...
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


% --- Executes just before UltrasoundImagingSimulator_Zheheng is made visible.
function UltrasoundImagingSimulator_Zheheng_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to UltrasoundImagingSimulator_Zheheng (see VARARGIN)

% Choose default command line output for UltrasoundImagingSimulator_Zheheng
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

mpth = mfilename('fullpath'); 
[mpth,~]=fileparts(mpth);
addpath(genpath(mpth));
field_init(0)

% global fs f0 c width element_height kerf focus N_elements image_width image_depth N_pixelx N_pixelz;
% ParametersSettoDefaultValues();

% UIWAIT makes UltrasoundImagingSimulator_Zheheng wait for user response (see UIRESUME)
% uiwait(handles.figure1);

function ParametersSettoDefaultValues()
nameanddefaultvalueoftheparameters=nameanddefaultvalueofparameters();
for indx=1:size(nameanddefaultvalueoftheparameters,1)
    nameanddefaultvalueoftheparameters{indx,2}=nameanddefaultvalueoftheparameters{indx,2}/unitconversiontoorderofmagnitude(nameanddefaultvalueoftheparameters{indx,1});
end
set(handles.parameterstable,'Data',nameanddefaultvalueoftheparameters);
% f0=5e6;                  %  Transducer center frequency [Hz]
% fs=100e6;                %  Sampling frequency [Hz]
% c=1540;                  %  Speed of sound [m/s]
% lambda=c/f0;             %  Wavelength [m]
% width=lambda/2;          %  Width of element
% element_height=7/1000;   %  Height of element [m]
% kerf=0.0025/1000;        %  Kerf [m]
% focus=[0 0 70];          %  Fixed focal point [m]
% N_elements=16;           %  Number of physical elements
% image_depth=0.005;       %  Depth of image along z axis
% image_width=0.005;       %  Width of image along x axis
% N_pixelx=128;            %  Number of pixels along x axis
% N_pixelz=128;            %  Number of pixels along z axis


function out1=nameanddefaultvalueofparameters() % Do not change the order of parameters
out1={'Sampling frequency (MHz)',50e6;...
    'Transducer center frequency (MHz)',5e6;...
    'Speed of sound (m/s)',1540;...
    'Width of element (mm)',1540/10e6;...
    'Height of element (mm)',0.007;...
    'Kerf (mm)',25e-7;...
    'Number of physical elements',32;...
    'Width of image along x axis (mm)',0.008;...
    'Depth of image along z axis (mm)',0.008;...
    'Number of pixels along x axis',64;...
    'Number of pixels along z axis',128;...
    'Focal length (mm)',10;
    'Distance from array to image (mm)',0.002};


% --- Outputs from this function are returned to the command line.
function varargout = UltrasoundImagingSimulator_Zheheng_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in exitbutton.
function exitbutton_Callback(hObject, eventdata, handles)
% hObject    handle to exitbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clear all
close all


function out1=unitconversiontoorderofmagnitude(in1)
if strcmp(in1(end-3:end),'(mm)')
    out1=1e-3;
elseif strcmp(in1(end-4:end),'(MHz)')
    out1=1e6;
else
    out1=1;
end

% --- Executes on selection change in parametermenu.
function parametermenu_Callback(hObject, eventdata, handles)
% hObject    handle to parametermenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns parametermenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from parametermenu
val_ParameterOptions = get(hObject,'Value');
nmndfltvlfprmtrs=get(handles.parameterstable,'Data');
set(handles.parametersedit,'string',num2str(nmndfltvlfprmtrs{val_ParameterOptions,2}))

% --- Executes during object creation, after setting all properties.
function parametermenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to parametermenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
nameanddefaultvalueoftheparameters=nameanddefaultvalueofparameters();
set(hObject,'String',nameanddefaultvalueoftheparameters(:,1));


function parametersedit_Callback(hObject, eventdata, handles)
% hObject    handle to parametersedit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of parametersedit as text
%        str2double(get(hObject,'String')) returns contents of parametersedit as a double


% --- Executes during object creation, after setting all properties.
function parametersedit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to parametersedit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject,'String','50');


% --- Executes on button press in comfirmsetbutton.
function comfirmsetbutton_Callback(hObject, eventdata, handles)
% hObject    handle to comfirmsetbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
val_parameter=str2double(get(handles.parametersedit,'string'));
if isnan(val_parameter)
    set(handles.parametersedit,'string','Error')
    return;
end
nvparameters=get(handles.parameterstable,'Data');
indx=get(handles.parametermenu,'Value');
nvparameters{indx,2}=val_parameter;
set(handles.parameterstable,'Data',nvparameters);
% switch get(handles.parametermenu,'Value')
%     case 1
%         fs=val_parameter*1e6;
%     case 2
%         f0=val_parameter*1e6;
%     case 3
%         c=val_parameter;
%     case 4
%         width=val_parameter/1e3;
%     case 5
%         element_height=val_parameter/1e3;
%     case 6
%         kerf=val_parameter/1e3;
%     case 7
%         N_elements=val_parameter;
%     case 8
%         image_width=val_parameter/1e3;
%     case 9
%         image_depth=val_parameter/1e3;
%     case 10
%         N_pixelx=val_parameter;
%     case 11
%         N_pixelz=val_parameter;
% end


% --- Executes on button press in Settodefaultbutton.
function Settodefaultbutton_Callback(hObject, eventdata, handles)
indx=get(handles.parametermenu,'Value');
nvpara=nameanddefaultvalueofparameters();
temp=nvpara{indx,2}/unitconversiontoorderofmagnitude(nvpara{indx,1});
nvpara=get(handles.parameterstable,'Data');
nvpara{indx,2}=temp;
set(handles.parametersedit,'string',num2str(nvpara{indx,2}));
set(handles.parameterstable,'Data',nvpara)
% switch get(handles.parametermenu,'Value')
%     case 1
%         fs=100e6;
%         set(handles.parametersedit,'string','100');
%     case 2
%         f0=5e6;
%         set(handles.parametersedit,'string','5');
%     case 3
%         c=1540;
%         set(handles.parametersedit,'string','1540');
%     case 4
%         width=1540/200e6;
%         set(handles.parametersedit,'string',num2str(width*1e3));
%     case 5
%         element_height=0.007;
%         set(handles.parametersedit,'string','7');
%     case 6
%         kerf=25e-7;
%         set(handles.parametersedit,'string','0.0025');
%     case 7
%         N_elements=16;
%         set(handles.parametersedit,'string','16');
%     case 8
%         image_width=0.005;
%         set(handles.parametersedit,'string','5');
%     case 9
%         image_depth=0.005;
%         set(handles.parametersedit,'string','5');
%     case 10
%         N_pixelx=128;
%         set(handles.parametersedit,'string','128');
%     case 11
%         N_pixelz=128;
%         set(handles.parametersedit,'string','128');
% end
% hObject    handle to Settodefaultbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function parameterstable_CreateFcn(hObject, eventdata, handles)
% hObject    handle to parameterstable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
nameanddefaultvalueoftheparameters=nameanddefaultvalueofparameters();
for indx=1:size(nameanddefaultvalueoftheparameters,1)
    nameanddefaultvalueoftheparameters{indx,2}=nameanddefaultvalueoftheparameters{indx,2}/unitconversiontoorderofmagnitude(nameanddefaultvalueoftheparameters{indx,1});
end
set(hObject,'Data',nameanddefaultvalueoftheparameters);


% --- Executes on button press in beamdisplaybutton.
function beamdisplaybutton_Callback(hObject, eventdata, handles)
% hObject    handle to beamdisplaybutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
nvpara=get(handles.parameterstable,'Data');
f0=nvpara{2,2}*unitconversiontoorderofmagnitude(nvpara{2,1});
fs=nvpara{1,2}*unitconversiontoorderofmagnitude(nvpara{1,1});
impulse_response=sin(2*pi*f0*(0:1/fs:2/f0));
impulse_response=impulse_response.*hanning(length(impulse_response))';
if get(handles.pulsecheckbox,'value')
    figure('name','Pulse')
    subplot(211)
    plot((0:1/fs:2/f0)*1e9,impulse_response);axis tight;grid on;xlabel('t/ns');title('Time response of the ultrasound pulse')
    subplot(212)
    plot(fs/3e8*(-150:150),fftshift(abs(fft(impulse_response,301))));axis tight;grid on;xlabel('f/MHz');title('Frequency response of the ultrasound pulse')
end
if ~get(handles.beamcheckbox,'value')
    return
end
set(hObject,'enable','off');
width=nvpara{4,2}*unitconversiontoorderofmagnitude(nvpara{4,1});
element_height=nvpara{5,2}*unitconversiontoorderofmagnitude(nvpara{5,1});
kerf=nvpara{6,2}*unitconversiontoorderofmagnitude(nvpara{6,1});
focus=[0 0 nvpara{12,2}*unitconversiontoorderofmagnitude(nvpara{12,1})];
N_elements=nvpara{7,2}*unitconversiontoorderofmagnitude(nvpara{7,1});
array_length=N_elements*(width+kerf)-kerf;
image_width=nvpara{8,2}*unitconversiontoorderofmagnitude(nvpara{8,1});
image_depth=nvpara{9,2}*unitconversiontoorderofmagnitude(nvpara{9,1});
dist_am=nvpara{13,2}*unitconversiontoorderofmagnitude(nvpara{13,1});
pixel_dx=image_width/round(256*abs(image_width)/max(abs([image_width image_depth])));
pixel_dz=image_depth/round(256*abs(image_depth)/max(abs([image_width image_depth])));
set_sampling(fs);
xmit_aperture = xdc_linear_array (N_elements, width, element_height, kerf, 1, 1, focus);
xdc_impulse(xmit_aperture, impulse_response);
% receive_aperture = xdc_linear_array (N_elements, width, element_height, kerf, 1, 1, focus);
[X_im,Z_im]=meshgrid((pixel_dx-image_width)/2:pixel_dx:image_width/2,dist_am+pixel_dz/2:pixel_dz:(image_depth+dist_am-pixel_dz/2));
[h,~]=calc_hp(xmit_aperture,[X_im(:),zeros(size(X_im(:))),Z_im(:)]);
h_indensity=sum(h.^2,1);
figure('name','Beam')
imagesc([-image_width/2 image_width/2]*1000,[dist_am image_depth+dist_am]*1000,10*log10(reshape(h_indensity,size(X_im))/max(h_indensity(:))),[lbdr(handles) 0]);
colorbar;colormap gray;axis equal;axis tight;
title('Normalised beam strength (dB)');xlabel('x(mm)');ylabel('z(mm)');
hold on
ax1=plot([-array_length/2 array_length/2]*1000,[0 0],'b','LineWidth',3);
legend(ax1,{'Transducer array'},'Location','bestoutside')
set(hObject,'enable','on');


function lbdisplayrangeedit_Callback(hObject, eventdata, handles)
% hObject    handle to lbdisplayrangeedit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lbdisplayrangeedit as text
%        str2double(get(hObject,'String')) returns contents of lbdisplayrangeedit as a double


% --- Executes during object creation, after setting all properties.
function lbdisplayrangeedit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lbdisplayrangeedit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject,'String','-40');

function out=lbdr(handles)
try
    out=str2double(get(handles.lbdisplayrangeedit,'string'));
catch
    set(handles.lbdisplayrangeedit,'string','Error');
    out=-40;
end
if isnan(out) || out>0 || out<-100
    set(handles.lbdisplayrangeedit,'string','Error');
    out=-40;
end


% --- Executes on button press in pulsecheckbox.
function pulsecheckbox_Callback(hObject, eventdata, handles)
% hObject    handle to pulsecheckbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of pulsecheckbox


% --- Executes on button press in beamcheckbox.
function beamcheckbox_Callback(hObject, eventdata, handles)
% hObject    handle to beamcheckbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of beamcheckbox


% --- Executes during object creation, after setting all properties.
function beamcheckbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to beamcheckbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
set(hObject,'value',1);



function phtfileedit_Callback(hObject, eventdata, handles)
% hObject    handle to phtfileedit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of phtfileedit as text
%        str2double(get(hObject,'String')) returns contents of phtfileedit as a double


% --- Executes during object creation, after setting all properties.
function phtfileedit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to phtfileedit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject,'string','8x_shrink_cyst_pht.mat')


% --- Executes on button press in displayphtpushbutton.
function displayphtpushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to displayphtpushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(hObject,'enable','off');
mpth = mfilename('fullpath'); 
[mpth,~]=fileparts(mpth);
phtfname=get(handles.phtfileedit,'string');
try
    load([mpth '/matfiles_ussim/' phtfname],'phantom_positions','phantom_amplitudes');
catch
    set(handles.phtfileedit,'string','Error');
    set(hObject,'enable','on');
    return;
end
% nvpara=get(handles.parameterstable,'Data');
% width=nvpara{4,2}*unitconversiontoorderofmagnitude(nvpara{4,1});
% phantom_positions(:,3)=phantom_positions(:,3)+width/2;
phantom_xlb=min(phantom_positions(:,1)); % min x
phantom_xub=max(phantom_positions(:,1)); % max x
phantom_zlb=min(phantom_positions(:,3)); % min z
phantom_zub=max(phantom_positions(:,3)); % max z
length_x=phantom_xub-phantom_xlb;
length_z=phantom_zub-phantom_zlb;
N_pixelx=round(512/max([length_z length_x])*length_x);
N_pixelz=round(512/max([length_z length_x])*length_z);
dx_im=length_x/N_pixelx;
dz_im=length_z/N_pixelz;
img0=zeros(N_pixelz,N_pixelx);
for indx=1:size(phantom_amplitudes,1)
    pp=[phantom_positions(indx,1)-phantom_xlb phantom_positions(indx,3)-phantom_zlb];
    i_xl=floor((pp(1)+dx_im/2)/dx_im);
    delxl=(pp(1)+dx_im/2)/dx_im-floor((pp(1)+dx_im/2)/dx_im);
    i_xu=i_xl+1;
    delxu=1-delxl;
    i_zl=floor((pp(2)+dz_im/2)/dz_im);
    delzl=(pp(2)+dz_im/2)/dz_im-floor((pp(2)+dz_im/2)/dz_im);
    i_zu=i_zl+1;
    delzu=1-delzl;
    if i_xl>0&&i_zl>0
        img0(i_zl,i_xl)=img0(i_zl,i_xl)+delxu*delzu*phantom_amplitudes(indx);
    end
    if i_zl>0&&i_xu<=N_pixelx
        img0(i_zl,i_xu)=img0(i_zl,i_xu)+delzu*delxl*phantom_amplitudes(indx);
    end
    if i_zu<=N_pixelz&&i_xl>0
        img0(i_zu,i_xl)=img0(i_zu,i_xl)+delzl*delxu*phantom_amplitudes(indx);
    end
    if i_xu<=N_pixelx&&i_zu<=N_pixelz
        img0(i_zu,i_xu)=img0(i_zu,i_xu)+delxl*delzl*phantom_amplitudes(indx);
    end
end
ampmax=max(img0(:));
ampmin=min(img0(:));
figure('name','Phantom')
imagesc([phantom_xlb phantom_xub]*1e3,[phantom_zlb phantom_zub]*1e3,img0);colormap(nonlincolormp_pht(ampmax,ampmin));
colorbar;axis tight equal;title('Phantom');xlabel('x(mm)');ylabel('z(mm)');
set(hObject,'enable','on');

function out=nonlincolormp_pht(amax,amin)
mmax=max(abs([amax amin]));
if amax<-0.01 || amin>0.01
    out=gray;
    return
end
mag_var=amin:(amax-amin)/200:amax;
nonlinmag_var=100/(mmax*log(10)*exp(2))*mag_var;
nonlinmag_var(abs(mag_var)>mmax*exp(2)/100)=sign(mag_var(abs(mag_var)>mmax*exp(2)/100)).*(log10(abs(mag_var(abs(mag_var)>mmax*exp(2)/100))/mmax)/2+1);
out=zeros(length(nonlinmag_var),3);
for indx=1:length(nonlinmag_var)
    if abs(nonlinmag_var(indx))<0.5
        out(indx,round([1.5 2.5]+0.5*sign(nonlinmag_var(indx))))=[1 1]-2*abs(nonlinmag_var(indx));
        out(indx,round(2-sign(nonlinmag_var(indx))))=1;
    else
        out(indx,round(2-sign(nonlinmag_var(indx))))=1-1.5*(abs(nonlinmag_var(indx))-0.5);
    end
end



function fblbedit_Callback(hObject, eventdata, handles)
% hObject    handle to fblbedit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fblbedit as text
%        str2double(get(hObject,'String')) returns contents of fblbedit as a double


% --- Executes during object creation, after setting all properties.
function fblbedit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fblbedit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject,'string','3.0294');



function fbubedit_Callback(hObject, eventdata, handles)
% hObject    handle to fbubedit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fbubedit as text
%        str2double(get(hObject,'String')) returns contents of fbubedit as a double


% --- Executes during object creation, after setting all properties.
function fbubedit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fbubedit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject,'string','7.0294')


% --- Executes on button press in genfbpushbutton.
function genfbpushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to genfbpushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
lcf=str2double(get(handles.fblbedit,'string'))*unitconversiontoorderofmagnitude(get(handles.fblbtext,'string'));
ucf=str2double(get(handles.fbubedit,'string'))*unitconversiontoorderofmagnitude(get(handles.fbubtext,'string'));
nvpara=get(handles.parameterstable,'Data');
f0=nvpara{2,2}*unitconversiontoorderofmagnitude(nvpara{2,1});
fs=nvpara{1,2}*unitconversiontoorderofmagnitude(nvpara{1,1});
bandwdth=abs(ucf-lcf);
impulse_response=sin(2*pi*f0*(0:1/fs:2/f0));
impulse_response=impulse_response.*hanning(length(impulse_response))';
impulse_spct=abs(fft(impulse_response,10001));
N_fre=floor(length(impulse_spct)/2)+1;
impulse_spct=impulse_spct(1:N_fre);
N_bf=round(N_fre*bandwdth*2/fs);
e_bd=conv(impulse_spct(:),ones([N_bf 1]));
e_bd=e_bd(N_bf:N_fre);
[~,lcf]=max(e_bd);
lcf=fs/2*lcf/N_fre;
set(handles.fblbedit,'string',num2str(lcf/unitconversiontoorderofmagnitude(get(handles.fblbtext,'string'))))
ucf=lcf+bandwdth;
set(handles.fbubedit,'string',num2str(ucf/unitconversiontoorderofmagnitude(get(handles.fbubtext,'string'))))


% --- Executes on button press in genrfpushbutton.
function genrfpushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to genrfpushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(hObject,'enable','off');
mpth = mfilename('fullpath'); 
[mpth,~]=fileparts(mpth);
phtfname=get(handles.phtfileedit,'string');
try
    load([mpth '/matfiles_ussim/' phtfname],'phantom_positions','phantom_amplitudes');
catch
    set(handles.phtfileedit,'string','Error');
    set(hObject,'enable','on');
    return;
end
nvpara=get(handles.parameterstable,'Data');
f0=nvpara{2,2}*unitconversiontoorderofmagnitude(nvpara{2,1});
fs=nvpara{1,2}*unitconversiontoorderofmagnitude(nvpara{1,1});
width=nvpara{4,2}*unitconversiontoorderofmagnitude(nvpara{4,1});
element_height=nvpara{5,2}*unitconversiontoorderofmagnitude(nvpara{5,1});
kerf=nvpara{6,2}*unitconversiontoorderofmagnitude(nvpara{6,1});
N_elements=nvpara{7,2}*unitconversiontoorderofmagnitude(nvpara{7,1});
focus=[0 0 nvpara{12,2}*unitconversiontoorderofmagnitude(nvpara{12,1})];
dist_am=nvpara{13,2}*unitconversiontoorderofmagnitude(nvpara{13,1});
set_sampling(fs);
xmit_aperture = xdc_linear_array (N_elements, width, element_height, kerf, 1, 1, focus);
impulse_response=sin(2*pi*f0*(0:1/fs:2/f0));
impulse_response=impulse_response.*hanning(length(impulse_response))';
xdc_impulse(xmit_aperture, impulse_response);
receive_aperture = xdc_linear_array (N_elements, width, element_height, kerf, 1, 1, focus);
phantom_positions(:,3)=phantom_positions(:,3)+dist_am;
[rf_data, tstart]=calc_scat_multi(xmit_aperture, receive_aperture, phantom_positions, phantom_amplitudes);
rf_data=single(rf_data);
rffilename=get(handles.rffilenameedit,'string');
try
    save([mpth '/matfiles_ussim/' rffilename],'rf_data','tstart');
catch
    set(handles.rffilenameedit,'string','Error, save in rf_data_temporary.mat');
    save([mpth '/matfiles_ussim/rf_data_temporary.mat'],'rf_data','tstart');
end
set(hObject,'enable','on');



function rffilenameedit_Callback(hObject, eventdata, handles)
% hObject    handle to rffilenameedit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rffilenameedit as text
%        str2double(get(hObject,'String')) returns contents of rffilenameedit as a double


% --- Executes during object creation, after setting all properties.
function rffilenameedit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rffilenameedit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject,'string','rf_data.mat');


% --- Executes on button press in genApushbutton.
function genApushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to genApushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(hObject,'enable','off');
mpth = mfilename('fullpath'); 
[mpth,~]=fileparts(mpth);
nvpara=get(handles.parameterstable,'Data');
f0=nvpara{2,2}*unitconversiontoorderofmagnitude(nvpara{2,1});%  Transducer center frequency [Hz]
fs=nvpara{1,2}*unitconversiontoorderofmagnitude(nvpara{1,1});%  Sampling frequency [Hz]
c=nvpara{3,2}*unitconversiontoorderofmagnitude(nvpara{3,1});%  Speed of sound [m/s]
width=nvpara{4,2}*unitconversiontoorderofmagnitude(nvpara{4,1});%  Width of element
element_height=nvpara{5,2}*unitconversiontoorderofmagnitude(nvpara{5,1});%  Height of element [m]
kerf=nvpara{6,2}*unitconversiontoorderofmagnitude(nvpara{6,1});%  Kerf [m]
N_elements=nvpara{7,2}*unitconversiontoorderofmagnitude(nvpara{7,1});%  Number of physical elements
focus=[0 0 nvpara{12,2}*unitconversiontoorderofmagnitude(nvpara{12,1})];%  Fixed focal point [m]
dist_am=nvpara{13,2}*unitconversiontoorderofmagnitude(nvpara{13,1});          
image_depth=nvpara{9,2}*unitconversiontoorderofmagnitude(nvpara{9,1});
image_width=nvpara{8,2}*unitconversiontoorderofmagnitude(nvpara{8,1});
N_pixelx=nvpara{10,2}*unitconversiontoorderofmagnitude(nvpara{10,1});
N_pixelz=nvpara{11,2}*unitconversiontoorderofmagnitude(nvpara{11,1});
pixel_dx=image_width/(N_pixelx-1);
pixel_dz=image_depth/(N_pixelz-1);
lcf=str2double(get(handles.fblbedit,'string'))*unitconversiontoorderofmagnitude(get(handles.fblbtext,'string'));
ucf=str2double(get(handles.fbubedit,'string'))*unitconversiontoorderofmagnitude(get(handles.fbubtext,'string'));
set_sampling(fs);
xmit_aperture = xdc_linear_array (N_elements, width, element_height, kerf, 1, 1, focus);
impulse_response=sin(2*pi*f0*(0:1/fs:2/f0));
impulse_response=impulse_response.*hanning(length(impulse_response))';
xdc_impulse(xmit_aperture, impulse_response);
receive_aperture = xdc_linear_array (N_elements, width, element_height, kerf, 1, 1, focus);
[X,Z]=meshgrid(-image_width/2:pixel_dx:image_width/2,dist_am:pixel_dz:image_depth+dist_am);
[rf_data_p,tstart_p]=calc_scat_multi(xmit_aperture, receive_aperture, [X(1) 0 Z(1)], 1);
rf_fre=rf_Ttof(rf_data_p,tstart_p,2*dist_am/c,2*(dist_am+image_depth)/c,lcf,ucf,fs);
A=single(zeros(2*length(rf_fre(:)),length(X(:))));
A(:,1)=single([real(rf_fre(:));imag(rf_fre(:))]);
for indx=2:length(X(:))
    [rf_data_p,tstart_p]=calc_scat_multi(xmit_aperture, receive_aperture, [X(indx) 0 Z(indx)], 1);
    rf_fre=rf_Ttof(rf_data_p,tstart_p,2*dist_am/c,2*(dist_am+image_depth)/c,lcf,ucf,fs);
    A(:,indx)=single([real(rf_fre(:));imag(rf_fre(:))]);
end
try
    save([mpth '/matfiles_ussim/' get(handles.genAedit,'string')],'A');
catch
    set(handles.genAedit,'string','Error, save in matrix_A_temporary.mat');
    save([mpth '/matfiles_ussim/matrix_A_temporary.mat'],'A');
end
set(hObject,'enable','on');


function out=rf_Ttof(rf,tstart,ts,te,lcf,ucf,fs)
if tstart>ts
    ns_temp=floor(fs*(tstart-ts));
    rf=conv2(rf,[zeros(ns_temp,1);1]);
    tstart=tstart-ns_temp/fs;
end
if tstart+size(rf,1)/fs<te
    rf=conv2(rf,[1;zeros(floor((te-tstart)*fs-size(rf,1)),1)]);
end
if tstart<ts
    ns_temp=ceil(fs*(ts-tstart));
    rf=rf((ns_temp+1):end,:);
    tstart=tstart+ns_temp/fs;
end
if tstart+size(rf,1)/fs>te
    rf=rf(1:floor(fs*(te-tstart)+1),:);
end
ns_temp=size(rf,1)-1;
rf=fft(rf,[],1);
if lcf>=0 && lcf<=ucf
    out=rf(round(lcf/fs*ns_temp+1):round(ucf/fs*ns_temp+1),:);
else
    out=rf(round(min(abs([lcf ucf]))/fs*ns_temp+1):round(max(abs([lcf ucf]))/fs*ns_temp+1),:);
end


function genAedit_Callback(hObject, eventdata, handles)
% hObject    handle to genAedit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of genAedit as text
%        str2double(get(hObject,'String')) returns contents of genAedit as a double


% --- Executes during object creation, after setting all properties.
function genAedit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to genAedit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject,'string','matrix_A.mat')


% --- Executes on button press in lsqrpushbutton.
function lsqrpushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to lsqrpushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(hObject,'enable','off');
mpth = mfilename('fullpath'); 
[mpth,~]=fileparts(mpth);
nvpara=get(handles.parameterstable,'Data');
fs=nvpara{1,2}*unitconversiontoorderofmagnitude(nvpara{1,1});%  Sampling frequency [Hz]
c=nvpara{3,2}*unitconversiontoorderofmagnitude(nvpara{3,1});%  Speed of sound [m/s]
dist_am=nvpara{13,2}*unitconversiontoorderofmagnitude(nvpara{13,1});          
image_depth=nvpara{9,2}*unitconversiontoorderofmagnitude(nvpara{9,1});
lcf=str2double(get(handles.fblbedit,'string'))*unitconversiontoorderofmagnitude(get(handles.fblbtext,'string'));
ucf=str2double(get(handles.fbubedit,'string'))*unitconversiontoorderofmagnitude(get(handles.fbubtext,'string'));
try
    load([mpth '/matfiles_ussim/' get(handles.rffilenameedit,'string')],'rf_data','tstart');
catch
    set(handles.rffilenameedit,'string','Error');
    set(hObject,'enable','on');
    return
end
noise_std=std_dev_noise(norm(rf_data(:))/sqrt(length(rf_data(:))),handles);
rf_data=rf_data+noise_std*randn(size(rf_data));
rf_data=rf_Ttof(rf_data,tstart,2*dist_am/c,2*(dist_am+image_depth)/c,lcf,ucf,fs);
rf_data=single([real(rf_data(:));imag(rf_data(:))]);
noise_std=std_dev_noise(double(norm(rf_data(:))/sqrt(length(rf_data(:)))),handles);
rf_data=rf_data+single(noise_std*randn(size(rf_data)));
try
    load([mpth '/matfiles_ussim/' get(handles.genAedit,'string')],'A');
catch
    set(handles.genAedit,'string','Error');
    set(hObject,'enable','on');
    return
end
n_iter_lsqr=str2double(get(handles.lsqrnumiteredit,'string'));
if isnan(n_iter_lsqr)
    n_iter_lsqr=8;
    set(handles.lsqrnumiteredit,'string','8');
end
img=lsqr_b(A,rf_data,n_iter_lsqr);
N_pixelx=nvpara{10,2}*unitconversiontoorderofmagnitude(nvpara{10,1});
N_pixelz=nvpara{11,2}*unitconversiontoorderofmagnitude(nvpara{11,1});
img=reshape(img(:,n_iter_lsqr),[N_pixelz N_pixelx]);
image_width=nvpara{8,2}*unitconversiontoorderofmagnitude(nvpara{8,1});
figure('name','LSQR reconstruction')
imagesc([-image_width/2 image_width/2]*1000,[dist_am image_depth+dist_am]*1000,20*log10(abs(img)/max(abs(img(:)))),[lbdr(handles) 0]);
colorbar;colormap gray;axis equal;axis tight;
title('LSQR reconstruction (dB)');xlabel('x(mm)');ylabel('z(mm)');
set(hObject,'enable','on');



function n=std_dev_noise(s,handles)
snr_rf=str2double(get(handles.snredit,'string'));
if isnan(snr_rf)
    snr_rf=10000;
    set(handles.snredit,'string','NaN');
end
n=s/10^(snr_rf/20);


function lsqrnumiteredit_Callback(hObject, eventdata, handles)
% hObject    handle to lsqrnumiteredit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lsqrnumiteredit as text
%        str2double(get(hObject,'String')) returns contents of lsqrnumiteredit as a double


% --- Executes during object creation, after setting all properties.
function lsqrnumiteredit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lsqrnumiteredit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject,'string','8');


% --- Executes during object creation, after setting all properties.
function lsqruipanel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lsqruipanel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in lsqrdefaultpushbutton.
function lsqrdefaultpushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to lsqrdefaultpushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.lsqrnumiteredit,'string','8');



function snredit_Callback(hObject, eventdata, handles)
% hObject    handle to snredit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of snredit as text
%        str2double(get(hObject,'String')) returns contents of snredit as a double


% --- Executes during object creation, after setting all properties.
function snredit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to snredit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject,'string','40')


% --- Executes on button press in yall1pushbutton.
function yall1pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to yall1pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(hObject,'enable','off');
mpth = mfilename('fullpath'); 
[mpth,~]=fileparts(mpth);
nvpara=get(handles.parameterstable,'Data');
fs=nvpara{1,2}*unitconversiontoorderofmagnitude(nvpara{1,1});%  Sampling frequency [Hz]
c=nvpara{3,2}*unitconversiontoorderofmagnitude(nvpara{3,1});%  Speed of sound [m/s]
dist_am=nvpara{13,2}*unitconversiontoorderofmagnitude(nvpara{13,1});          
image_depth=nvpara{9,2}*unitconversiontoorderofmagnitude(nvpara{9,1});
lcf=str2double(get(handles.fblbedit,'string'))*unitconversiontoorderofmagnitude(get(handles.fblbtext,'string'));
ucf=str2double(get(handles.fbubedit,'string'))*unitconversiontoorderofmagnitude(get(handles.fbubtext,'string'));
try
    load([mpth '/matfiles_ussim/' get(handles.rffilenameedit,'string')],'rf_data','tstart');
catch
    set(handles.rffilenameedit,'string','Error');
    set(hObject,'enable','on');
    return
end
noise_std=std_dev_noise(norm(rf_data(:))/sqrt(length(rf_data(:))),handles);
rf_data=rf_data+noise_std*randn(size(rf_data));
rf_data=rf_Ttof(rf_data,tstart,2*dist_am/c,2*(dist_am+image_depth)/c,lcf,ucf,fs);
rf_data=single([real(rf_data(:));imag(rf_data(:))]);
noise_std=std_dev_noise(norm(double(rf_data(:))/sqrt(length(rf_data(:)))),handles);
rf_data=rf_data+single(noise_std*randn(size(rf_data)));
try
    load([mpth '/matfiles_ussim/' get(handles.genAedit,'string')],'A');
catch
    set(handles.genAedit,'string','Error');
    set(hObject,'enable','on');
    return
end
opts.nonorth=1;
opts.rho=str2double(get(handles.rhoedit,'string'));
if isnan(opts.rho)
    set(handles.rhoedit,'string','Error');
    opts.rho=0.001;
end
opts.tol=str2double(get(handles.yall1toledit,'string'));
if isnan(opts.tol)
    set(handles.yall1toledit,'string','Error');
    opts.tol=0.01;
end
N_pixelx=nvpara{10,2}*unitconversiontoorderofmagnitude(nvpara{10,1});
N_pixelz=nvpara{11,2}*unitconversiontoorderofmagnitude(nvpara{11,1});
dwtmode('per');
if get(handles.l1regubs1radiobutton,'value')
    flag_bs=[get(handles.idmatcheckbox,'value');get(handles.dct2checkbox,'value');...
        get(handles.Haarcheckbox,'value');get(handles.meycheckbox,'value');...
        get(handles.meyxHaarzcheckbox,'value');get(handles.Haarxmeyzbox,'value')];
    opts.basis.times=@(x) transform_func(x,flag_bs,1,N_pixelx,N_pixelz);
    opts.basis.trans=@(x) transform_func(x,flag_bs,0,N_pixelx,N_pixelz);
end
if get(handles.l1regubs2radiobutton,'value')
    filenm=get(handles.bshandletxtedit,'string');
    for ind=1:length(filenm)
        if (filenm(end)<'a'||filenm(end)>'z')&&(filenm(end)<'A'||filenm(end)>'Z')
            filenm(end)=[];
        else
            break;
        end
    end
    if length(filenm)>4 && strcmp(filenm(end-3:end),'.mat')
        load([mpth '/matfiles_ussim/' filenm],'phi');
        opts.basis.times=@(x) phi*x;
        opts.basis.trans=@(x) phi'*x;
    else
        handletxt=importdata([mpth '/' filenm]);
        eval(['opts.basis.times=' handletxt{1} ';']);
        eval(['opts.basis.trans=' handletxt{2} ';']);
    end
end
opts.print=2;
img=yall1(1e20*A,1e20*rf_data,opts);
img=reshape(img(:),[N_pixelz N_pixelx]);
image_width=nvpara{8,2}*unitconversiontoorderofmagnitude(nvpara{8,1});
figure('name','YALL1 reconstruction')
imagesc([-image_width/2 image_width/2]*1000,[dist_am image_depth+dist_am]*1000,20*log10(abs(img)/max(abs(img(:)))),[lbdr(handles) 0]);
colorbar;colormap gray;axis equal;axis tight;
title('YALL1 reconstruction (dB)');xlabel('x(mm)');ylabel('z(mm)');
set(hObject,'enable','on');

function out=transform_func(x,flag_bs,flag_notinv,nx,nz)
if sum(flag_bs(2:end))<0.5
    out=x;
    return;
end
if flag_notinv
out=[];
if flag_bs(1)
    out=x;
end
if flag_bs(2)
    temp=reshape(x,[nz nx]);
    temp=dct2(temp);
    out=[out;temp(:)];
end
fn={'haar','dmey'};
for indx=3:4
    if flag_bs(indx)
        temp=reshape(x,[nz nx]);
        temp=conv2(temp,[1;zeros(mod(nz,2),1)]);
        temp=conv2(temp,[1 zeros(1,mod(nx,2))]);
        [temp,temp1,temp2,temp3]=dwt2(temp,fn{indx-2});
        temp=[temp,temp2;temp1,temp3];
        out=[out;temp(:)];
    end
end
for indx=5:6
    if flag_bs(indx)
        temp=reshape(x,[nz nx]);
        temp=conv2(temp,[1;zeros(mod(nz,2),1)]);
        temp=conv2(temp,[1 zeros(1,mod(nx,2))]);
        for ind=1:size(temp,1)
            [temp1,temp2]=dwt(temp(ind,:),fn{7-indx});
            temp(ind,:)=[temp1(:)',temp2(:)'];
        end
        for ind=1:size(temp,2)
            [temp1,temp2]=dwt(temp(:,ind),fn{indx-4});
            temp(:,ind)=[temp1(:);temp2(:)];
        end
        out=[out;temp(:)];
    end
end
out=out/sqrt(sum(flag_bs));
else
out=0;    
if flag_bs(1)
    out=x(1:nz*nx);
    x(1:nz*nx)=[];
end
if flag_bs(2)
    temp=reshape(x(1:nz*nx),[nz nx]);
    x(1:nz*nx)=[];
    temp=idct2(temp);
    out=out+temp(:);
end
fn={'haar','dmey'};
for indx=3:4
    if flag_bs(indx)
        temp1=nz+mod(nz,2);
        temp2=nx+mod(nx,2);
        temp=reshape(x(1:temp1*temp2),[temp1 temp2]);
        temp=idwt2(temp(1:temp1/2,1:temp2/2),temp(temp1/2+1:end,1:temp2/2),temp(1:temp1/2,temp2/2+1:end),temp(temp1/2+1:end,temp2/2+1:end),fn{indx-2});
        temp=temp(1:nz,1:nx);
        out=out+temp(:);
    end
end
for indx=5:6
    if flag_bs(indx)
        temp1=nz+mod(nz,2);
        temp2=nx+mod(nx,2);
        temp=reshape(x(1:temp1*temp2),[temp1 temp2]);
        for ind=1:size(temp,2)
            temp3=idwt(temp(1:temp1/2,ind),temp(temp1/2+1:end,ind),fn{indx-4});
            temp(:,ind)=temp3(:);
        end
        for ind=1:size(temp,1)
            temp3=idwt(temp(ind,1:temp2/2),temp(ind,temp2/2+1:end),fn{7-indx});
            temp(ind,:)=temp3(:)';
        end
        temp=temp(1:nz,1:nx);
        out=out+temp(:);
    end
end
out=out/sqrt(sum(flag_bs));
end

function rhoedit_Callback(hObject, eventdata, handles)
% hObject    handle to rhoedit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rhoedit as text
%        str2double(get(hObject,'String')) returns contents of rhoedit as a double


% --- Executes during object creation, after setting all properties.
function rhoedit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rhoedit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject,'string','0.001');



function yall1toledit_Callback(hObject, eventdata, handles)
% hObject    handle to yall1toledit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of yall1toledit as text
%        str2double(get(hObject,'String')) returns contents of yall1toledit as a double


% --- Executes during object creation, after setting all properties.
function yall1toledit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to yall1toledit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject,'string','1e-2');


% --- Executes on button press in l1regubs1radiobutton.
function l1regubs1radiobutton_Callback(hObject, eventdata, handles)
% hObject    handle to l1regubs1radiobutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of l1regubs1radiobutton


function bshandletxtedit_Callback(hObject, eventdata, handles)
% hObject    handle to bshandletxtedit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bshandletxtedit as text
%        str2double(get(hObject,'String')) returns contents of bshandletxtedit as a double


% --- Executes during object creation, after setting all properties.
function bshandletxtedit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bshandletxtedit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject,'string','transform_handles.m');


% --- Executes on button press in dct2checkbox.
function dct2checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to dct2checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of dct2checkbox



% --- Executes on button press in l1regubs2radiobutton.
function l1regubs2radiobutton_Callback(hObject, eventdata, handles)
% hObject    handle to l1regubs2radiobutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of l1regubs2radiobutton


% --- Executes on button press in idmatcheckbox.
function idmatcheckbox_Callback(hObject, eventdata, handles)
% hObject    handle to idmatcheckbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of idmatcheckbox


% --- Executes on button press in Haarcheckbox.
function Haarcheckbox_Callback(hObject, eventdata, handles)
% hObject    handle to Haarcheckbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Haarcheckbox


% --- Executes on button press in meycheckbox.
function meycheckbox_Callback(hObject, eventdata, handles)
% hObject    handle to meycheckbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of meycheckbox


% --- Executes on button press in Haarxmeyzbox.
function Haarxmeyzbox_Callback(hObject, eventdata, handles)
% hObject    handle to Haarxmeyzbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Haarxmeyzbox


% --- Executes on button press in meyxHaarzcheckbox.
function meyxHaarzcheckbox_Callback(hObject, eventdata, handles)
% hObject    handle to meyxHaarzcheckbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of meyxHaarzcheckbox
