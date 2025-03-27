function varargout = OrbitSimulator(varargin)
% ORBITSIMULATOR MATLAB code for OrbitSimulator.fig
%      ORBITSIMULATOR, by itself, creates a new ORBITSIMULATOR or raises the existing
%      singleton
%
%      H = ORBITSIMULATOR returns the handle to a new ORBITSIMULATOR or the handle to
%      the existing singleton*.
%
%      ORBITSIMULATOR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ORBITSIMULATOR.M with the given input arguments.
%
%      ORBITSIMULATOR('Property','Value',...) creates a new ORBITSIMULATOR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before OrbitSimulator_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to OrbitSimulator_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools meT_anomaly.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
% Edit the above text to modify the response to help OrbitSimulator
% Last Modified by GUIDE v2.5 06-Jan-2017 14:43:47
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @OrbitSimulator_OpeningFcn, ...
                   'gui_OutputFcn',  @OrbitSimulator_OutputFcn, ...
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
% --- Executes just before OrbitSimulator is made visible.
function OrbitSimulator_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to OrbitSimulator (see VARARGIN)
% Choose default command line output for OrbitSimulator
handles.output = hObject;
% Update handles structure
guidata(hObject, handles);
% UIWAIT makes OrbitSimulator wait for user response (see UIRESUME)
% uiwait(handles.figure1);
% --- Outputs from this function are returned to the command line.
function varargout = OrbitSimulator_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Get default command line output from handles structure
varargout{1} = handles.output;
% --- Executes on button press in browse.
function browse_Callback(hObject, eventdata, handles)
[filename, pathname] = uigetfile('*.txt');
set( handles.direct, 'String', [pathname filename],'Enable','on');
% hObject    handle to browse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- Executes on selection change in method.
function method_Callback(hObject, eventdata, handles)
% hObject    handle to method (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: contents = cellstr(get(hObject,'String')) returns method contents as cell array
%        contents{get(hObject,'Value')} returns selected item from method
% --- Executes during object creation, after setting all properties.
function method_CreateFcn(hObject, eventdata, handles)
% hObject    handle to method (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --- Executes on key press with focus on browse and none of its controls.
function browse_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to browse (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
function direct_Callback(hObject, eventdata, handles)
% hObject    handle to direct (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of direct as text
%        str2double(get(hObject,'String')) returns contents of direct as a double
% --- Executes during object creation, after setting all properties.
function direct_CreateFcn(hObject, eventdata, handles)
% hObject    handle to direct (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
switch get(handles.popupmenu2,'Value')
   
    case 2
        set(handles.rx,'Visible','on')
        set(handles.ry,'Visible','on')
        set(handles.rz,'Visible','on')
        set(handles.vx,'Visible','on')
        set(handles.vy,'Visible','on')
        set(handles.vz,'Visible','on')
        set(handles.in1,'Visible','on')
        set(handles.in2,'Visible','on')
        set(handles.in3,'Visible','on')
        set(handles.in4,'Visible','on')
        set(handles.in5,'Visible','on')
        set(handles.in6,'Visible','on')
        set(handles.in7,'Visible','off')
        set(handles.in8,'Visible','off')
        set(handles.in9,'Visible','off')
        set(handles.in10,'Visible','off')
        set(handles.in11,'Visible','off')
        set(handles.in12,'Visible','off')
        set(handles.browse,'Visible','off')
        set(handles.direct,'Visible','off')
        set(handles.a,'Visible','off')
        set(handles.i,'Visible','off')
        set(handles.raan,'Visible','off')
        set(handles.e,'Visible','off')
        set(handles.aop,'Visible','off')
        set(handles.ta,'Visible','off')
    
    case 3
        set(handles.a,'Visible','on')
        set(handles.i,'Visible','on')
        set(handles.raan,'Visible','on')
        set(handles.e,'Visible','on')
        set(handles.aop,'Visible','on')
        set(handles.ta,'Visible','on')
        set(handles.rx,'Visible','off')
        set(handles.ry,'Visible','off')
        set(handles.rz,'Visible','off')
        set(handles.vx,'Visible','off')
        set(handles.vy,'Visible','off')
        set(handles.vz,'Visible','off')
        set(handles.in1,'Visible','off')
        set(handles.in2,'Visible','off')
        set(handles.in3,'Visible','off')
        set(handles.in4,'Visible','off')
        set(handles.in5,'Visible','off')
        set(handles.in6,'Visible','off')
        set(handles.in7,'Visible','on')
        set(handles.in8,'Visible','on')
        set(handles.in9,'Visible','on')
        set(handles.in10,'Visible','on')
        set(handles.in11,'Visible','on')
        set(handles.in12,'Visible','on')
        set(handles.browse,'Visible','off')
        set(handles.direct,'Visible','off')
    case 4
        set(handles.direct,'Visible','on')
        set(handles.browse,'Visible','on')
        set(handles.rx,'Visible','off')
        set(handles.ry,'Visible','off')
        set(handles.rz,'Visible','off')
        set(handles.vx,'Visible','off')
        set(handles.vy,'Visible','off')
        set(handles.vz,'Visible','off')
        set(handles.in1,'Visible','off')
        set(handles.in2,'Visible','off')
        set(handles.in3,'Visible','off')
        set(handles.in4,'Visible','off')
        set(handles.in5,'Visible','off')
        set(handles.in6,'Visible','off')
        set(handles.in7,'Visible','off')
        set(handles.in8,'Visible','off')
        set(handles.in9,'Visible','off')
        set(handles.in10,'Visible','off')
        set(handles.in11,'Visible','off')
        set(handles.in12,'Visible','off')
        set(handles.a,'Visible','off')
        set(handles.i,'Visible','off')
        set(handles.raan,'Visible','off')
        set(handles.e,'Visible','off')
        set(handles.aop,'Visible','off')
        set(handles.ta,'Visible','off')
        
    otherwise
        set(handles.direct,'Visible','off')
        set(handles.browse,'Visible','off')
        set(handles.rx,'Visible','off')
        set(handles.ry,'Visible','off')
        set(handles.rz,'Visible','off')
        set(handles.vx,'Visible','off')
        set(handles.vy,'Visible','off')
        set(handles.vz,'Visible','off')
        set(handles.in1,'Visible','off')
        set(handles.in2,'Visible','off')
        set(handles.in3,'Visible','off')
        set(handles.in4,'Visible','off')
        set(handles.in5,'Visible','off')
        set(handles.in6,'Visible','off')
        set(handles.in7,'Visible','off')
        set(handles.in8,'Visible','off')
        set(handles.in9,'Visible','off')
        set(handles.in10,'Visible','off')
        set(handles.in11,'Visible','off')
        set(handles.in12,'Visible','off')
        set(handles.a,'Visible','off')
        set(handles.i,'Visible','off')
        set(handles.raan,'Visible','off')
        set(handles.e,'Visible','off')
        set(handles.aop,'Visible','off')
        set(handles.ta,'Visible','off')
end
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2
% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function in1_Callback(hObject, eventdata, handles)
get(handles.in1,'string');
% hObject    handle to in1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of in1 as text
%        str2double(get(hObject,'String')) returns contents of in1 as a double
% --- Executes during object creation, after setting all properties.
function in1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to in1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function in2_Callback(hObject, eventdata, handles)
get(handles.in2,'string');
% hObject    handle to in2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of in2 as text
%        str2double(get(hObject,'String')) returns contents of in2 as a double
% --- Executes during object creation, after setting all properties.
function in2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to in2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function in3_Callback(hObject, eventdata, handles)
get(handles.in3,'string');
% hObject    handle to in3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of in3 as text
%        str2double(get(hObject,'String')) returns contents of in3 as a double
% --- Executes during object creation, after setting all properties.
function in3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to in3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function in4_Callback(hObject, eventdata, handles)
get(handles.in4,'string');
% hObject    handle to in4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of in4 as text
%        str2double(get(hObject,'String')) returns contents of in4 as a double
function in5_Callback(hObject, eventdata, handles)
get(handles.in5,'string');
% hObject    handle to in5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of in5 as text
%        str2double(get(hObject,'String')) returns contents of in5 as a double
function in6_Callback(hObject, eventdata, handles)
get(handles.in6,'string');
% hObject    handle to in6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of in6 as text
%        str2double(get(hObject,'String')) returns contents of in6 as a double
% --- Executes on button press in calc.
function calc_Callback(hObject, eventdata, handles)
%clearing the axis for new calculation%
cla(handles.axes5)
cla(handles.axes4)
 rev_no=str2double(get(handles.rev,'String')); %extracting the number of revolutions
 earth_rad=6371000;
 if isnan(rev_no)==1        %error for null input of the revolutions around earth
    errordlg('You must specify the revolutions number','error');
 end
 switch get(handles.method,'Value')
     case 1
         errordlg('You must choose a calculation method','error');
         switch get(handles.popupmenu2,'Value')
                case 1
                    errordlg('You must choose an input method','error')
         end
        case 2
        mu = 398.6004418e12; 
        earth_rad=6371000;
        runspeed=str2double(get(handles.runspeed,'String'));
        if str2double(get(handles.runspeed,'String'))==0
            runspeed=300;
        end
        
        switch get(handles.popupmenu2,'Value')
            case 1
                errordlg('You must choose an input method','error')
            case 2
                mu = 398.6004418e12; 
                rx=str2double(get(handles.in1,'String'));
                ry=str2double(get(handles.in2,'String'));
                rz=str2double(get(handles.in3,'String'));
                vx=str2double(get(handles.in4,'String'));
                vy=str2double(get(handles.in5,'String'));
                vz=str2double(get(handles.in6,'String'));
                r=[rx ry rz];
                v=[vx vy vz];
                v=v*1000;
                r=r*1000;
                rmag = sqrt(dot(r, r));
                vmag=sqrt(dot(v, v));
                E=(vmag^2/2)-(mu/rmag);
                a=(-mu/(2*E));
                hv=cross(r,v);
                hmag = sqrt(hv(1,1)^2+hv(1,2)^2+hv(1,3)^2);
                hhat = hv/hmag;
                p=hmag^2/mu;
                e=sqrt(1-(p/a));
                D=dot(r, v);
            if D>0
                T_anomaly=acosd((p/(e*rmag))-(1/e));
            else
                T_anomaly=360-(acosd((p/(e*rmag))-(1/e)));
            end
            inclination=acosd(hv(1,3)/hmag);
            k=[0 0 1];
            i=[1 0 0];
            n_vec=cross(k,hv);
            n=sqrt(n_vec(1,1)^2+n_vec(1,2)^2+n_vec(1,3)^2);
            if n_vec(2)>0
                RAANd=acosd(dot(n_vec,i)/n);
            else
                RAANd=360-(acosd(dot(n_vec,i)/n));
            end
                d=dot(r,v);
                f=vmag^2-(mu/rmag);
                ecc=(1/mu)*((f*r)-(d*v));
            if ecc(3)>0
                AOP=acosd((dot(ecc,n_vec))/(e*n));
            else 
                AOP=360-(acosd((dot(ecc,n_vec))/(e*n)));
        end
            case 3
                T_anomaly=str2double(get(handles.in10,'String'));
                RAANd=str2double(get(handles.in8,'String'));
                inclination=str2double(get(handles.in7,'String'));
                AOP=str2double(get(handles.in9,'String'));
                a=str2double(get(handles.in11,'String'))*1000;
                e=str2double(get(handles.in12,'String'));
                p=a*(1-e^2);
                rmag=(a*(1-e^2))/(1+(e*cosd(T_anomaly)));
                vmag=sqrt(mu*((2/rmag)-(1/a)));
                rPQW=rmag*[cosd(T_anomaly) sind(T_anomaly) 0];
                vPQW=sqrt(mu/p)*[-sind(T_anomaly) e+cosd(T_anomaly) 0];
                Rijk=[(cosd(RAANd)*cosd(AOP))-(sind(RAANd)*sind(AOP)*cosd(inclination)) -(cosd(RAANd)*sind(AOP))-(sind(RAANd)*cosd(AOP)*cosd(inclination)) sind(RAANd)*sind(inclination);(sind(RAANd)*cosd(AOP))+(cosd(RAANd)*sind(AOP)*cosd(inclination)) -(sind(RAANd)*sind(AOP))+(cosd(RAANd)*cosd(AOP)*cosd(inclination)) -cosd(RAANd)*sind(inclination);(sind(AOP)*sind(inclination)) (cosd(AOP)*sind(inclination)) cosd(inclination)];
                Vijk=[(cosd(RAANd)*cosd(AOP))-(sind(RAANd)*sind(AOP)*cosd(inclination)) -(cosd(RAANd)*sind(AOP))-(sind(RAANd)*cosd(AOP)*cosd(inclination)) sind(RAANd)*sind(inclination);(sind(RAANd)*cosd(AOP))+(cosd(RAANd)*sind(AOP)*cosd(inclination)) -(sind(RAANd)*sind(AOP))+(cosd(RAANd)*cosd(AOP)*cosd(inclination)) -cosd(RAANd)*sind(inclination);(sind(AOP)*sind(inclination)) (cosd(AOP)*sind(inclination)) cosd(inclination)];
                r= Rijk*transpose(rPQW);
                v= Vijk*transpose(vPQW);
                r=r';
                v=v';
            case 4
            % TLE file name 
            fname = get(handles.direct,'String');
            % Open the TLE file and read TLE elements
            fid = fopen(fname, 'rb');
            while 1
            % read second line
            tline = fgetl(fid);
            if ~ischar(tline), break, end
            inclination = str2double(tline(9:16));                                 % Orbit Inclination (degrees)
            RAANd= str2double(tline(18:25));                               % Right Ascension of Ascending Node (degrees)
            e = str2double(strcat('0.',tline(27:33)));                   % Eccentricity
            AOP = str2double(tline(35:42));                                  % Argument of Perigee (degrees)
            M = str2double(tline(44:51));                                  % Mean Anomaly (degrees) 
            a = ( (398600/(str2double(tline(53:63))*2*pi/86400)^2 )^(1/3))*1000 ;   % semi major axis 
            end
            fclose(fid);
            Md= M*(pi/180);
            f=@(E)(Md+e*sin(E)-E);
            EA=fsolve(f,0);
            EA=double(EA);
            T_anomaly=acosd((cos(EA)-e)/(1-(e*cos(EA))));
                p=a*(1-e^2);
                rmag=(a*(1-e^2))/(1+(e*cosd(T_anomaly)));
                vmag=sqrt(mu*((2/rmag)-(1/a)));
                rPQW=rmag*[cosd(T_anomaly) sind(T_anomaly) 0];
                vPQW=sqrt(mu/p)*[-sind(T_anomaly) e+cosd(T_anomaly) 0];
                Rijk=[(cosd(RAANd)*cosd(AOP))-(sind(RAANd)*sind(AOP)*cosd(inclination)) -(cosd(RAANd)*sind(AOP))-(sind(RAANd)*cosd(AOP)*cosd(inclination)) sind(RAANd)*sind(inclination);(sind(RAANd)*cosd(AOP))+(cosd(RAANd)*sind(AOP)*cosd(inclination)) -(sind(RAANd)*sind(AOP))+(cosd(RAANd)*cosd(AOP)*cosd(inclination)) -cosd(RAANd)*sind(inclination);(sind(AOP)*sind(inclination)) (cosd(AOP)*sind(inclination)) cosd(inclination)];
                Vijk=[(cosd(RAANd)*cosd(AOP))-(sind(RAANd)*sind(AOP)*cosd(inclination)) -(cosd(RAANd)*sind(AOP))-(sind(RAANd)*cosd(AOP)*cosd(inclination)) sind(RAANd)*sind(inclination);(sind(RAANd)*cosd(AOP))+(cosd(RAANd)*sind(AOP)*cosd(inclination)) -(sind(RAANd)*sind(AOP))+(cosd(RAANd)*cosd(AOP)*cosd(inclination)) -cosd(RAANd)*sind(inclination);(sind(AOP)*sind(inclination)) (cosd(AOP)*sind(inclination)) cosd(inclination)];
                r= Rijk*transpose(rPQW);
                v= Vijk*transpose(vPQW);
                r=r';
                v=v';
        end
        J2000_days=103752/24;
        rhat = r/rmag; 
        vhat = v/vmag;
        V=v;
        %angles in radians
        v= T_anomaly*(pi/180);
        w= AOP*(pi/180);
        RAAN= RAANd*(pi/180);
        inc= inclination*(pi/180);
        %MEAN ANOMALY (M0)
        M0 = 2*atan(sqrt((1-e)/(1+e))*tan(v/2)) - e*sqrt(1-e^2)*sin(v)/(1+e*cos(v)); %in rad
        %ERROR if Launch position is inside the Earth
        if (sqrt (r(1)*r(1)+r(2)*r(2)+r(3)*r(3)) <= 6731000)
            blast (r(1), r(2), r(3), 2000000);
            errordlg('Launch Position Inside Earth','Launch error');
            return;
        end
        %Final Adjustments to Initial Mean Anomaly
        if M0<0
            M0=-M0;
        end
            E=M0;
            for i=1:5
                E = E + (M0 + e*sin(E) - E)/(1 - e*cos(E));
            end
            v= 2*atan(sqrt((1+e)/(1-e))*tan(E/2));
            R = a*(1-e*cos(E));
            Xeci = R*(cos(w + v)*cos(RAAN) - sin(w+v)*sin(RAAN)*cos(inc));
            Yeci = R*(cos(w + v)*sin(RAAN) + sin(w+v)*cos(RAAN)*cos(inc));
            Zeci = R*(sin(w + v)*sin(inc));
            c=0;
            while abs(r(1)-Xeci)>100 && abs(r(2)-Yeci)>100 && abs(r(3)-Zeci)>100 && c<15
                if c~=0
                    if (c<3)
                        M0=2*pi-M0;
                    end
                    if (c<5 && c>=3)
                        M0=-M0;
                    end
                    if (c<10 && c>=5)
                        M0=M0+(pi/2);
                        if (c==9)
                            M0=iniM0;
                        end
                    end
                    if (c>=10 && c<15)
                        M0=M0-(pi/2);
                        if (c==15)
                            M0=iniM0;
                        end
                    end 
                end
                E=M0;
                for i=1:5
                    E = E + (M0 + e*sin(E) - E)/(1 - e*cos(E));
                end
                v= 2*atan(sqrt((1+e)/(1-e))*tan(E/2));
                R = a*(1-e*cos(E));
                b  = a*sqrt(1-e^2);
                Xeci = R*(cos(w + v)*cos(RAAN) - sin(w+v)*sin(RAAN)*cos(inc));
                Yeci = R*(cos(w + v)*sin(RAAN) + sin(w+v)*cos(RAAN)*cos(inc));
                Zeci = R*(sin(w + v)*sin(inc));
                c=c+1;
        end
        if (M0>2*pi)
            M0=M0-2*pi;
        end
        if (M0<0)
            M0=M0+2*pi;     
        end
        %Orbital Period 
        orbital_period = sqrt((a*a*a*4*pi*pi)/mu);
        %ERROR if eccentricity is greater than supported
        if (e<0)
            if (e>0.95)
                errordlg('Eccentricity greater than supported','eccentricity error')
            end
        end
            if (e>1)
                errordlg('Eccentricity greater than supported (e>1)','eccentricity error')
            return
            end
            
        check_error=[r v];
        check=zeros(1,length(check_error));
        if isequal(isnan(check_error),check)==0;
            errordlg('Your inputs are missing please recheck them','error')
        else
        set(handles.axes5,'Visible','on')
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %DRAWING THE STATIC VISUALIZATION COMPONENTS%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Initializing the Drawing Space
        axes(handles.axes5)
        lim=(1+e)*a; %Setting the limits of the graph
        axis([-lim, lim, -lim, lim, -lim, lim])	
        view(150,15) 
        axis equal
        if get(handles.scalebox,'Value')==1
            axis on
        else
            axis off
        end
        shg
        hold on
        rotate3d on
        if get(handles.gridbox,'Value')==1
            grid on
        end
points=[lim,lim,-lim,-lim;lim,-lim,-lim,lim;0,0,0,0];
equatorial=fill3(points(1,:),points(2,:),points(3,:),'b');
orbitalp=fill3(points(1,:),points(2,:),points(3,:),'c');
rotate(orbitalp,[0,1,0],RAANd)
rotate(orbitalp,[1,0,0],inclination)
alpha(0.1)
if get(handles.equapbox,'Value')==1
  set(equatorial,'Visible','on');
else
    set(equatorial,'Visible','off');
end
if get(handles.orbitpbox,'Value')==1
  set(orbitalp,'Visible','on');
else
    set(orbitalp,'Visible','off');
end
        %Plotting the Earth
        equat_rad=6378137.00;
            polar_rad=6356752.3142;
            [xx yy zz]=ellipsoid (0,0,0,equat_rad, equat_rad, polar_rad);
            load('topo.mat','topo','topomap1');
            topo2 = [topo(:,181:360) topo(:,1:180)];
            pro.FaceColor= 'texture';
            pro.EdgeColor = 'none';
            pro.FaceLighting = 'phong';
            pro.Cdata = topo2;
            earth= surface(xx,yy,zz,pro);
            colormap(topomap1)
        omega_earth = 7.292115855377074e-005; % (rad/sec)  
        Go = 1.727564365843028; % (rad)  
        GMST = Go + omega_earth*86400*(J2000_days + 0.5);
        GMST = GMST - 2*pi*floor(GMST/(2*pi));
        GMST_deg=GMST*(180/pi);
        rotate (earth, [0 0 1], GMST_deg);
         Xaxis= line([0 lim],[0 0],[0 0],'Color', 'red', 'Marker','.','LineStyle','-','Visible','on');
        Yaxis= line([0 0],[0 lim],[0 0],'Color', 'red', 'Marker','.','LineStyle','-','Visible','on');
        rotate (Xaxis, [0 0 1], GMST_deg);
        rotate (Yaxis, [0 0 1], GMST_deg);
        if get(handles.axisbox,'Value')==1
            set(Xaxis,'Visible','on');
            set(Yaxis,'Visible','on');
        else
            set(Xaxis,'Visible','off');
            set(Yaxis,'Visible','off');
        end
        %Plotting the ECI Axes
        ECIx=line([0 lim],[0 0],[0 0],'Color', 'black', 'Marker','.','LineStyle','-');
        ECIy=line([0 0],[0 lim],[0 0],'Color', 'black', 'Marker','.','LineStyle','-');
        ECIz=line([0 0],[0 0],[0 lim],'Color', 'black', 'Marker','.','LineStyle','-');
        if get(handles.axis2box,'Value')==1
            set(ECIx,'Visible','on')
            set(ECIy,'Visible','on')
            set(ECIz,'Visible','on')
        else
            set(ECIx,'Visible','off')
            set(ECIy,'Visible','off')
            set(ECIz,'Visible','off')
        end
        %Plotting Initial Velocity Vector
        vin=line([r(1) r(1)+2000*V(1)],[r(2) r(2)+2000*V(2)],[r(3) r(3)+2000*V(3)],'Color', 'green','Marker','.','LineWidth', 2, 'MarkerSize', 8,'LineStyle','-');
        %Plotting the initial poisition of the satellite
       pin=plot3 (r(1), r(2), r(3),'o', 'MarkerEdgeColor', 'black','MarkerFaceColor','green','MarkerSize', 10);
if get(handles.vinbox,'Value')==1
    set(vin,'Visible','on')
else
    set(vin,'Visible','off')
end
if get(handles.pinbox,'Value')==1
    set(pin,'Visible','on')
else
    set(pin,'Visible','off')
end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %DRAWING THE DYNAMIC VISUALIZATION COMPONENTS%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        k=0;
        long=1:rev_no*ceil(orbital_period/runspeed);
        %Plotting the movement of the satellite
        Xcoord(1)=r(1);
        Ycoord(1)=r(2);
        Zcoord(1)=r(3);
                        set(handles.axes4,'Visible','on');
        axes(handles.axes4);
        hold on
        image([0 360],[90 -90],topo,'CDataMapping', 'scaled');
        colormap(topomap1);
        axis equal
        axis ([0 360 -90 90]);
        cc=1;
        colorcount=1;
        ground=[];
        orbitline=[];
        for time = 1:rev_no*ceil(orbital_period/runspeed)
            k=k+1;
            %Computing Eccentric Anomaly
            axes(handles.axes5)
            if get(handles.gridbox,'Value')==1
                grid on
            else
                grid off
            end
        if get(handles.axisbox,'Value')==1
            set(Xaxis,'Visible','on');
            set(Yaxis,'Visible','on');
        else
            set(Xaxis,'Visible','off');
            set(Yaxis,'Visible','off');
        end
        if get(handles.axis2box,'Value')==1
            set(ECIx,'Visible','on')
            set(ECIy,'Visible','on')
            set(ECIz,'Visible','on')
        else
            set(ECIx,'Visible','off')
            set(ECIy,'Visible','off')
            set(ECIz,'Visible','off')
        end
        if get(handles.equapbox,'Value')==1
  set(equatorial,'Visible','on');
else
    set(equatorial,'Visible','off');
        end
if get(handles.orbitpbox,'Value')==1
  set(orbitalp,'Visible','on');
else
    set(orbitalp,'Visible','off');
end
        if get(handles.vinbox,'Value')==1
    set(vin,'Visible','on')
else
    set(vin,'Visible','off')
end
if get(handles.pinbox,'Value')==1
    set(pin,'Visible','on')
else
    set(pin,'Visible','off')
end
if get(handles.scalebox,'Value')==1
            axis on
        else
            axis off
        end
                    
            E=M0;
            for i=1:5
                E = E + (M0 + e*sin(E) - E)/(1 - e*cos(E));
            end
            %Computing the True Anomaly
            v= 2*atan(sqrt((1+e)/(1-e))*tan(E/2));
            %Computing 'r' in polar coordinates
            r = a*(1-e*cos(E));
            GMST=GMST+(runspeed/86400)*2*pi;
            %Computes the Cartesian Co-ordinates in ECI frame from 'r' and orbital
            %elements
            Xeci = r*(cos(w + v)*cos(RAAN) - sin(w+v)*sin(RAAN)*cos(inc));
            Yeci = r*(cos(w + v)*sin(RAAN) + sin(w+v)*cos(RAAN)*cos(inc));
            Zeci = r*(sin(w + v)*sin(inc));
            rotate (earth, [0 0 1], (runspeed)*(360/86400))
            rotate (Xaxis, [0 0 1], (runspeed)*(360/86400))
            rotate (Yaxis, [0 0 1], (runspeed)*(360/86400))
            %Drawing the red sphere
            array(k)=plot3 (Xeci, Yeci, Zeci,'o', 'MarkerEdgeColor', 'k','MarkerFaceColor','r','MarkerSize', 6);
            position(k)=line([0 Xeci],[0 Yeci], [0 Zeci],'Color', 'yellow', 'LineWidth', 2);
            if (k~=1)
            set (array(k-1), 'Visible', 'off');
            set (position(k-1), 'Visible', 'off');
            end
                if get(handles.posbox,'Value')==0
                    set (position(k), 'Visible', 'off');
                end
            if (time~=1 && time<=ceil(orbital_period/runspeed)+1)
                Xcoord(k)=Xeci;
                Ycoord(k)=Yeci;
                Zcoord(k)=Zeci;
                orbitline(time)=line([Xcoord(k-1) Xcoord(k)],[Ycoord(k-1) Ycoord(k)], [Zcoord(k-1) Zcoord(k)],'Color', 'black', 'LineWidth', 2);
            end
            if get(handles.orbitbox,'Value')==1 && time~=1
                set(orbitline(:),'Visible','on')
            else
                set(orbitline(:),'Visible','off')
            end
            if (GMST>2*pi)
                GMST=GMST-2*pi;
            end
            lat(k)=atan(Zeci/sqrt(Xeci*Xeci+Yeci*Yeci))*(180/pi);
            ECIX=[cos(GMST) sin(GMST) 0];
            Pos=[Xeci Yeci 0];
            cvec = cross(ECIX,Pos);
            angleyz = mod(sign(dot([0 0 1],cvec))*atan2(norm(cvec),dot(ECIX,Pos)),2*pi);
            long(k) =(180/pi)* angleyz;
            
                         pause (0.01);
            %Blast condition
            if (sqrt (Xeci*Xeci+Yeci*Yeci+Zeci*Zeci) <= 6731000)
                blast (Xeci, Yeci, Zeci, 2000000);
                errordlg('Satellite crashed into earth','Orbit intersects with earth');
                break;
            end
            if time==ceil(orbital_period/runspeed)*cc;
                colorcount=colorcount+1;
                cc=cc+1;
            end
            color=['r' 'c' 'm' 'y' 'w' 'b' 'k' 'g' 'r'];
            axes(handles.axes4)
            if (time~=1 && abs(long(k-1)-long(k))<100)
                ground(time)=line([long(k-1) long(k)],[-lat(k-1) -lat(k)],'Color', color(colorcount), 'LineWidth', 2.5);
            end
            
            h=plot (long(k),-lat(k),'o', 'MarkerEdgeColor', 'k','MarkerFaceColor','r','MarkerSize', 8);
            if get(handles.groundbox,'Value')==1 && time>1 && abs(long(k-1)-long(k))<100;
               set(ground(:),'Visible','on')
           else
               set(ground(:),'Visible','off')
           end
            pause(0.0001)
            if time~=rev_no*ceil(orbital_period/runspeed);
                delete(h);
            end
            M0=M0+sqrt(mu/(a*a*a))*runspeed; %Updating Mean Anomaly for next iteration
            if colorcount==9;
                colorcount=1;
            end
            vPQW=sqrt(398.6004418e12/p)*[-sin(v) e+cos(v) 0];
                Vijk=[(cosd(RAANd)*cosd(AOP))-(sind(RAANd)*sind(AOP)*cosd(inclination)) -(cosd(RAANd)*sind(AOP))-(sind(RAANd)*cosd(AOP)*cosd(inclination)) sind(RAANd)*sind(inclination);(sind(RAANd)*cosd(AOP))+(cosd(RAANd)*sind(AOP)*cosd(inclination)) -(sind(RAANd)*sind(AOP))+(cosd(RAANd)*cosd(AOP)*cosd(inclination)) -cosd(RAANd)*sind(inclination);(sind(AOP)*sind(inclination)) (cosd(AOP)*sind(inclination)) cosd(inclination)];
                V= Vijk*transpose(vPQW);
                v1=num2str(V(1,1)/1000);
                v2=num2str(V(2,1)/1000);
                v3=num2str(V(3,1)/1000);
               r1=num2str(floor(Xeci/1000));
               r2=num2str(floor(Yeci/1000));
               r3=num2str(floor(Zeci/1000));
                set(handles.rxf,'String',['Rx= ' r1 ' km'])
                set(handles.ryf,'String',['Ry= ' r2 ' km'])
                set(handles.rzf,'String',['Rz= ' r3 ' km'])
                set(handles.vxf,'String',['Vx= ' v1 ' km/s'])
                set(handles.vyf,'String',['Vy= ' v2 ' km/s'])
                set(handles.vzf,'String',['Vz= ' v3 ' km/s'])
        end
 end
        case 3
            switch get(handles.popupmenu2,'Value')
                case 1
                    errordlg('You must choose an input method','error')
                case 2 
                t=500;
                rx=str2double(get(handles.in1,'String'));
                ry=str2double(get(handles.in2,'String'));
                rz=str2double(get(handles.in3,'String'));
                vx=str2double(get(handles.in4,'String'));
                vy=str2double(get(handles.in5,'String'));
                vz=str2double(get(handles.in6,'String'));
                rv=[rx ry rz vx vy vz];
                case 3
                    T_anomaly=str2double(get(handles.in10,'String'));
                    RAANd=str2double(get(handles.in8,'String'));
                    inclination=str2double(get(handles.in7,'String'));
                    AOP=str2double(get(handles.in9,'String'));
                    a=str2double(get(handles.in11,'String'));
                    e=str2double(get(handles.in12,'String'));
                    mu=398600;
                    p=a*(1-e^2);
                rmag=(a*(1-e^2))/(1+(e*cosd(T_anomaly)));
                vmag=sqrt(mu*((2/rmag)-(1/a)));
                rPQW=rmag*[cosd(T_anomaly) sind(T_anomaly) 0];
                vPQW=sqrt(mu/p)*[-sind(T_anomaly) e+cosd(T_anomaly) 0];
                Rijk=[(cosd(RAANd)*cosd(AOP))-(sind(RAANd)*sind(AOP)*cosd(inclination)) -(cosd(RAANd)*sind(AOP))-(sind(RAANd)*cosd(AOP)*cosd(inclination)) sind(RAANd)*sind(inclination);(sind(RAANd)*cosd(AOP))+(cosd(RAANd)*sind(AOP)*cosd(inclination)) -(sind(RAANd)*sind(AOP))+(cosd(RAANd)*cosd(AOP)*cosd(inclination)) -cosd(RAANd)*sind(inclination);(sind(AOP)*sind(inclination)) (cosd(AOP)*sind(inclination)) cosd(inclination)];
                Vijk=[(cosd(RAANd)*cosd(AOP))-(sind(RAANd)*sind(AOP)*cosd(inclination)) -(cosd(RAANd)*sind(AOP))-(sind(RAANd)*cosd(AOP)*cosd(inclination)) sind(RAANd)*sind(inclination);(sind(RAANd)*cosd(AOP))+(cosd(RAANd)*sind(AOP)*cosd(inclination)) -(sind(RAANd)*sind(AOP))+(cosd(RAANd)*cosd(AOP)*cosd(inclination)) -cosd(RAANd)*sind(inclination);(sind(AOP)*sind(inclination)) (cosd(AOP)*sind(inclination)) cosd(inclination)];
                r= Rijk*transpose(rPQW);
                v= Vijk*transpose(vPQW);
                r=r';
                v=v';
                rv=[r v];
                case 4
            % TLE file name 
            fname = get(handles.direct,'string');
            % Open the TLE file and read TLE elements
            fid = fopen(fname, 'rb');
            while 1
            % read second line
            tline = fgetl(fid);
            if ~ischar(tline), break, end
            inclination = str2double(tline(9:16));                                 % Orbit Inclination (degrees)
            RAANd= str2double(tline(18:25));                               % Right Ascension of Ascending Node (degrees)
            e = str2double(strcat('0.',tline(27:33)));                   % Eccentricity
            AOP = str2double(tline(35:42));                                  % Argument of Perigee (degrees)
            M = str2double(tline(44:51));                                  % Mean Anomaly (degrees) 
            a = ( (398600/(str2double(tline(53:63))*2*pi/86400)^2 )^(1/3))*1000 ;   % semi major axis 
            end
            fclose(fid);
            Md= M*(pi/180);
f=@(E)(Md+e*sin(E)-E);
            EA=fsolve(f,0);
            EA=double(EA);
            T_anomaly=acosd((cos(EA)-e)/(1-(e*cos(EA))));
                p=a*(1-e^2);
                rmag=(a*(1-e^2))/(1+(e*cosd(T_anomaly)));
                mu=398600*10^9;
                vmag=sqrt(mu*((2/rmag)-(1/a)));
                rPQW=rmag*[cosd(T_anomaly) sind(T_anomaly) 0];
                vPQW=sqrt(mu/p)*[-sind(T_anomaly) e+cosd(T_anomaly) 0];
                Rijk=[(cosd(RAANd)*cosd(AOP))-(sind(RAANd)*sind(AOP)*cosd(inclination)) -(cosd(RAANd)*sind(AOP))-(sind(RAANd)*cosd(AOP)*cosd(inclination)) sind(RAANd)*sind(inclination);(sind(RAANd)*cosd(AOP))+(cosd(RAANd)*sind(AOP)*cosd(inclination)) -(sind(RAANd)*sind(AOP))+(cosd(RAANd)*cosd(AOP)*cosd(inclination)) -cosd(RAANd)*sind(inclination);(sind(AOP)*sind(inclination)) (cosd(AOP)*sind(inclination)) cosd(inclination)];
                Vijk=[(cosd(RAANd)*cosd(AOP))-(sind(RAANd)*sind(AOP)*cosd(inclination)) -(cosd(RAANd)*sind(AOP))-(sind(RAANd)*cosd(AOP)*cosd(inclination)) sind(RAANd)*sind(inclination);(sind(RAANd)*cosd(AOP))+(cosd(RAANd)*sind(AOP)*cosd(inclination)) -(sind(RAANd)*sind(AOP))+(cosd(RAANd)*cosd(AOP)*cosd(inclination)) -cosd(RAANd)*sind(inclination);(sind(AOP)*sind(inclination)) (cosd(AOP)*sind(inclination)) cosd(inclination)];
                r= Rijk*transpose(rPQW);
                v= Vijk*transpose(vPQW);
                r=r'/1000;
                v=v'/1000;
                rv=[r v];
            end
            check_error=rv;
        check=zeros(1,length(check_error));
        if isequal(isnan(check_error),check)==0;
            errordlg('Your inputs are missing please recheck them','error')
        else
            t=0;
            orbit(t,rv);
                axes(handles.axes5);
                 initial=rv;
runspeed=str2double(get(handles.runspeed,'String'));
if str2double(get(handles.runspeed,'String'))==0
    runspeed=300;
end
r=[rv(1) rv(2) rv(3)]*1000;
v=[rv(4) rv(5) rv(6)]*1000;
[e,a,inclination,RAANd,AOP,T_anomaly,Tp]=orbitelem(r,v);
tend=rev_no*Tp;
RAANd=360-RAANd
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
J2000_days=103752/24;
V=v;
%ERROR if eccentricity is greater than supported
if (e<0)
    if (e>0.95)
        errordlg('Eccentricity greater than supported','eccentricity error')
    end
end
    if (e>1)
        errordlg('Eccentricity greater than supported (e>1)','eccentricity error')
        return
    end
   tspan=[0:runspeed:tend];
   [T,RV]=ode45(@orbit,tspan,initial);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DRAWING THE STATIC VISUALIZATION COMPONENTS%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Initializing the Drawing Space
      
set(handles.axes5,'Visible','on')
        axes(handles.axes5)
        lim=(1+e)*a; %Setting the limits of the graph
        axis([-lim, lim, -lim, lim, -lim, lim])	
        view(150,15) 
        axis equal
        shg
        rotate3d on
        if get(handles.scalebox,'Value')==1
                axis on
            else
                axis off
            end
        hold on
        grid on
         if get(handles.gridbox,'Value')==1
                grid on
            else
                grid off
            end
points=[lim,lim,-lim,-lim;lim,-lim,-lim,lim;0,0,0,0];
equatorial=fill3(points(1,:),points(2,:),points(3,:),'b');
orbitalp=fill3(points(1,:),points(2,:),points(3,:),'c');
rotate(orbitalp,[0,1,0],RAANd)
rotate(orbitalp,[1,0,0],inclination)
alpha(0.1)
if get(handles.equapbox,'Value')==1
  set(equatorial,'Visible','on');
else
    set(equatorial,'Visible','off');
end
if get(handles.orbitpbox,'Value')==1
  set(orbitalp,'Visible','on');
else
    set(orbitalp,'Visible','off');
end
%Plotting the Earth
equat_rad=6378137.00;
    polar_rad=6356752.3142;
    [xx yy zz]=ellipsoid (0,0,0,equat_rad, equat_rad, polar_rad);
    load('topo.mat','topo','topomap1');
    topo2 = [topo(:,181:360) topo(:,1:180)];
    pro.FaceColor= 'texture';
    pro.EdgeColor = 'none';
    pro.FaceLighting = 'phong';
    pro.Cdata = topo2;
    earth= surface(xx,yy,zz,pro);
    colormap(topomap1)
omega_earth = 7.292115855377074e-005; % (rad/sec)  
Go = 1.727564365843028; % (rad)  
GMST = Go + omega_earth*86400*(J2000_days + 0.5);
GMST = GMST - 2*pi*floor(GMST/(2*pi));
GMST_deg=GMST*(180/pi);
rotate (earth, [0 0 1], GMST_deg);
Xaxis= line([0 lim],[0 0],[0 0],'Color', 'red', 'Marker','.','LineStyle','-');
Yaxis= line([0 0],[0 lim],[0 0],'Color', 'red', 'Marker','.','LineStyle','-');
rotate (Xaxis, [0 0 1], GMST_deg);
rotate (Yaxis, [0 0 1], GMST_deg);
%Plotting the ECI Axes
ECIx=line([0 lim],[0 0],[0 0],'Color', 'black', 'Marker','.','LineStyle','-');
ECIy=line([0 0],[0 lim],[0 0],'Color', 'black', 'Marker','.','LineStyle','-');
ECIz=line([0 0],[0 0],[0 lim],'Color', 'black', 'Marker','.','LineStyle','-');
%Plotting Initial Velocity Vector
vin=line([r(1) r(1)+2000*V(1)],[r(2) r(2)+2000*V(2)],[r(3) r(3)+2000*V(3)],'Color', 'green','Marker','.','LineWidth', 2, 'MarkerSize', 8,'LineStyle','-');
%Plotting the initial poisition of the satellite
pin=plot3 (r(1), r(2), r(3),'o', 'MarkerEdgeColor', 'black','MarkerFaceColor','green','MarkerSize', 10);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DRAWING THE DYNAMIC VISUALIZATION COMPONENTS%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  set(handles.axes4,'Visible','on');
        axes(handles.axes4);
        hold on
        image([0 360],[90 -90],topo,'CDataMapping', 'scaled');
        colormap(topomap1);
        axis equal
        axis ([0 360 -90 90]);
%Plotting the movement of the satellite
Xcoord(1)=r(1);
Ycoord(1)=r(2);
Zcoord(1)=r(3);
k=0;
[j,~]=size(RV);
j=floor(j/rev_no);
cc=1;
colorcount=1;
           
        if get(handles.axisbox,'Value')==1
            set(Xaxis,'Visible','on');
            set(Yaxis,'Visible','on');
        else
            set(Xaxis,'Visible','off');
            set(Yaxis,'Visible','off');
        end
        if get(handles.axis2box,'Value')==1
            set(ECIx,'Visible','on')
            set(ECIy,'Visible','on')
            set(ECIz,'Visible','on')
        else
            set(ECIx,'Visible','off')
            set(ECIy,'Visible','off')
            set(ECIz,'Visible','off')
        end
        if get(handles.vinbox,'Value')==1
            set(vin,'Visible','on')
        else
            set(vin,'Visible','off')
        end
if get(handles.pinbox,'Value')==1
    set(pin,'Visible','on')
else
    set(pin,'Visible','off')
end
ground=[];
orbitline=[];
for i=2:length(RV(:,1))
    k=k+1;
    axes(handles.axes5)
            if get(handles.gridbox,'Value')==1
                grid on
            else
                grid off
            end
            if get(handles.scalebox,'Value')==1
                axis on
            else
                axis off
            end
            
        if get(handles.axisbox,'Value')==1
            set(Xaxis,'Visible','on');
            set(Yaxis,'Visible','on');
        else
            set(Xaxis,'Visible','off');
            set(Yaxis,'Visible','off');
        end
        if get(handles.axis2box,'Value')==1
            set(ECIx,'Visible','on')
            set(ECIy,'Visible','on')
            set(ECIz,'Visible','on')
        else
            set(ECIx,'Visible','off')
            set(ECIy,'Visible','off')
            set(ECIz,'Visible','off')
        end
        if get(handles.equapbox,'Value')==1
  set(equatorial,'Visible','on');
else
    set(equatorial,'Visible','off');
        end
if get(handles.orbitpbox,'Value')==1
  set(orbitalp,'Visible','on');
else
    set(orbitalp,'Visible','off');
end
        if get(handles.vinbox,'Value')==1
            set(vin,'Visible','on')
        else
            set(vin,'Visible','off')
        end
if get(handles.pinbox,'Value')==1
    set(pin,'Visible','on')
else
    set(pin,'Visible','off')
end
    GMST=GMST+(runspeed/86400)*2*pi;
    rotate (earth, [0 0 1], (runspeed)*(360/86400))
    rotate (Xaxis, [0 0 1], (runspeed)*(360/86400))
    rotate (Yaxis, [0 0 1], (runspeed)*(360/86400))
    Xeci=RV(i,1)*1000;
    Yeci=RV(i,2)*1000;
    Zeci=RV(i,3)*1000;
    array(k)=plot3 (Xeci, Yeci, Zeci,'o', 'MarkerEdgeColor', 'k','MarkerFaceColor','r','MarkerSize', 6);
    position(k)=line([0 Xeci],[0 Yeci], [0 Zeci],'Color', 'yellow', 'LineWidth', 2);
    if (k~=1)
    set (array(k-1), 'Visible', 'off');
    set (position(k-1), 'Visible', 'off');
    end
    if get(handles.posbox,'Value')==0
        set (position(k), 'Visible', 'off');
    end
if i~=2;
         Xcoord(k)=Xeci;
        Ycoord(k)=Yeci;
        Zcoord(k)=Zeci;
        orbitline(i)=line([Xcoord(k-1) Xcoord(k)],[Ycoord(k-1) Ycoord(k)], [Zcoord(k-1) Zcoord(k)],'Color', 'black', 'LineWidth', 2);
end
            if get(handles.orbitbox,'Value')==1 && i~=2
                set(orbitline(:),'Visible','on')
            else
                set(orbitline(:),'Visible','off')
            end
         if (GMST>2*pi)
        GMST=GMST-2*pi;
         end
         if i==j*cc;
             colorcount=colorcount+1;
             cc=cc+1;
         end
    color=['r' 'c' 'm' 'y' 'w' 'b' 'k' 'g' 'r'];
    
    lat(k)=atan(Zeci/sqrt(Xeci*Xeci+Yeci*Yeci))*(180/pi);
    ECIX=[cos(GMST) sin(GMST) 0];
    Pos=[Xeci Yeci 0];
    cvec = cross(ECIX,Pos);
    angleyz = mod(sign(dot([0 0 1],cvec))*atan2(norm(cvec),dot(ECIX,Pos)),2*pi);
    long(k) =(180/pi)* angleyz;
        if (sqrt (Xeci*Xeci+Yeci*Yeci+Zeci*Zeci) <= 6731000)
                blast (Xeci, Yeci, Zeci, 2000000);
                errordlg('Satellite crashed into earth','Orbit intersects with earth');
                break;
        end
        axes(handles.axes4)
            if (i~=2 && abs(long(k-1)-long(k))<100)
                ground(i)=line([long(k-1) long(k)],[-lat(k-1) -lat(k)],'Color',color(colorcount), 'LineWidth', 2.5);
            end
            h=plot (long(k),-lat(k),'o', 'MarkerEdgeColor', 'k','MarkerFaceColor','r','MarkerSize', 8);
           if get(handles.groundbox,'Value')==1 && i>2 && abs(long(k-1)-long(k))<100;
               set(ground(:),'Visible','on')
           else
               set(ground(:),'Visible','off')
           end
            pause(0.0001)
            if i~=length(RV(:,1))
                delete(h);
            end
            if colorcount==9;
                colorcount=1;
            end
                set(handles.rxf,'String',['Rx= ' num2str(floor(RV(i,1))) ' km'])
                set(handles.ryf,'String',['Ry= ' num2str(floor(RV(i,2))) ' km'])
                set(handles.rzf,'String',['Rz= ' num2str(floor(RV(i,3))) ' km'])
                set(handles.vxf,'String',['Vx= ' num2str(RV(i,4)) ' km/s'])
                set(handles.vyf,'String',['Vy= ' num2str(RV(i,5)) ' km/s'])
                set(handles.vzf,'String',['Vz= ' num2str(RV(i,6)) ' km/s'])
end
        end
 end
% hObject    handle to calc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
function rev_Callback(hObject, eventdata, handles)
get(handles.rev,'string');
% hObject    handle to rev (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of rev as text
%        str2double(get(hObject,'String')) returns contents of rev as a double
% --- Executes during object creation, after setting all properties.
function rev_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rev (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --- Executes during object creation, after setting all properties.
function in4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to in4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --- Executes during object creation, after setting all properties.
function in5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to in5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --- Executes during object creation, after setting all properties.
function in6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to in6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% hObject    handle to direct (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
function in7_Callback(hObject, eventdata, handles)
% hObject    handle to in7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of in7 as text
%        str2double(get(hObject,'String')) returns contents of in7 as a double
% --- Executes during object creation, after setting all properties.
function in7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to in7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function in8_Callback(hObject, eventdata, handles)
% hObject    handle to in8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of in8 as text
%        str2double(get(hObject,'String')) returns contents of in8 as a double
% --- Executes during object creation, after setting all properties.
function in8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to in8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function in9_Callback(hObject, eventdata, handles)
% hObject    handle to in9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of in9 as text
%        str2double(get(hObject,'String')) returns contents of in9 as a double
% --- Executes during object creation, after setting all properties.
function in9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to in9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function in10_Callback(hObject, eventdata, handles)
% hObject    handle to in10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of in10 as text
%        str2double(get(hObject,'String')) returns contents of in10 as a double
% --- Executes during object creation, after setting all properties.
function in10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to in10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function in11_Callback(hObject, eventdata, handles)
% hObject    handle to in11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of in11 as text
%        str2double(get(hObject,'String')) returns contents of in11 as a double
% --- Executes during object creation, after setting all properties.
function in11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to in11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function in12_Callback(hObject, eventdata, handles)
% hObject    handle to in12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of in12 as text
%        str2double(get(hObject,'String')) returns contents of in12 as a double
% --- Executes during object creation, after setting all properties.
function in12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to in12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit19_Callback(hObject, eventdata, handles)
% hObject    handle to rev (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of rev as text
%        str2double(get(hObject,'String')) returns contents of rev as a double
% --- Executes during object creation, after setting all properties.
function edit19_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rev (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% --- Executes during object creation, after setting all properties.
function axes6_CreateFcn(hObject, eventdata, handles)
I=imread('t.jpg');
imagesc(I);
colormap gray;
set(gca,'Visible','off');
% hObject    handle to axes6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: place code in OpeningFcn to populate axes6
% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over direct.
function direct_ButtonDownFcn(hObject, eventdata, handles)
set(hObject,'String','','Enable','on');
uicontrol(handles.direct);
% hObject    handle to direct (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- Executes during object creation, after setting all properties.
function axes7_CreateFcn(hObject, eventdata, handles)
I=imread('zewail.jpg');
imagesc(I);
set(gca,'Visible','off');
% hObject    handle to axes7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: place code in OpeningFcn to populate axes7
% --- Executes during object creation, after setting all properties.
function axes5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: place code in OpeningFcn to populate axes5
% --- Executes during object creation, after setting all properties.
function axes4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: place code in OpeningFcn to populate axes4
function runspeed_Callback(hObject, eventdata, handles)
% hObject    handle to runspeed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of runspeed as text
%        str2double(get(hObject,'String')) returns contents of runspeed as a double
% --- Executes during object creation, after setting all properties.
function runspeed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to runspeed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
set(handles.runspeed,'String','300')
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- Executes on button press in pushbutton5.
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over runspeed.
function runspeed_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to runspeed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% hObject    handle to runspeed (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
count=num2str(str2double(get(handles.runspeed,'String'))+50);
set(handles.runspeed,'String',count);
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
count=str2double(get(handles.runspeed,'String'))-50;
if count<=0
    count=0;
end
count=num2str(count);
set(handles.runspeed,'String',count);
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- Executes on button press in gridbox.
function gridbox_Callback(hObject, eventdata, handles)
% hObject    handle to gridbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of gridbox
% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over gridbox.
function gridbox_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to gridbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- Executes on button press in axisbox.
function axisbox_Callback(hObject, eventdata, handles)
% hObject    handle to axisbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of axisbox
% --- Executes on button press in axis2box.
function axis2box_Callback(hObject, eventdata, handles)
% hObject    handle to axis2box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of axis2box
% --- Executes on button press in vinbox.
function vinbox_Callback(hObject, eventdata, handles)
% hObject    handle to vinbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of vinbox
% --- Executes on button press in pinbox.
function pinbox_Callback(hObject, eventdata, handles)
% hObject    handle to pinbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of pinbox
% --- Executes on button press in posbox.
function posbox_Callback(hObject, eventdata, handles)
% hObject    handle to posbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of posbox
% --- Executes on button press in scalebox.
function scalebox_Callback(hObject, eventdata, handles)
if get(handles.scalebox,'Value')==1
    set(handles.gridbox,'Enable','on')
else
    set(handles.gridbox,'Enable','off')
end
% hObject    handle to scalebox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of scalebox
% --- Executes on button press in groundbox.
function groundbox_Callback(hObject, eventdata, handles)
% hObject    handle to groundbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of groundbox
% --- Executes on button press in orbitbox.
function orbitbox_Callback(hObject, eventdata, handles)
% hObject    handle to orbitbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of orbitbox
% --- Executes on button press in equapbox.
function equapbox_Callback(hObject, eventdata, handles)
% hObject    handle to equapbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of equapbox
% --- Executes on button press in orbitpbox.
function orbitpbox_Callback(hObject, eventdata, handles)
% hObject    handle to orbitpbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of orbitpbox