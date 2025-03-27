function [POLAR, TARGET] = mission_init()
% add all subfolders to path
try
    Xi([0;0;0;1]);
catch
    addpath("Algorithms/")
    addpath("Subsystems/")
    addpath("common_utils/data")
    addpath("common_utils/data/igrf")
    addpath("common_utils/time/")                
    addpath("common_utils/frames/")
    addpath("common_utils/math/")
end

% Dynamic Parameters
Ixx = .052; Iyy = .05; Izz = .025;
Ixy = 1e-4; Ixz = 1e-4; Iyz = 1e-4;
POLAR.I = [Ixx Ixy Ixz; Ixy Iyy Iyz; Ixz Iyz Izz ]; % Inertia Matrix of spacecraft [kg m^2]
POLAR.m = 10;                                       % mass of POLAR [kg]
POLAR.C_d = 2.2;                                    % drag coefficient


% Set Initial Conditions
a = 6780;  %  altitude [km]
e = 1e-4;  %  eccentricity [unitless]
i = 51;    %  inclination  [degrees]
RAAN = 1.9;  %  RAAN         [degrees]
argp = -.6;  %  Arg. of perigee [degrees]
M = .1;    %  Mean Anomaly [radians]
POLAR.OE = [a,e,i,RAAN,argp,M];                   % orbital elements of POLAR 
[r,v] = kep2rv(a,e,i,argp,RAAN,M);                                                                
POLAR.r = r; POLAR.v = v;                         % initiate position and velocity of POLAR in ECI frame
q = [1 1 2 8]'; POLAR.q = q/norm(q);          % set random quaternion, and normalize it
POLAR.omega = [0.2;0.4;0.3]*(pi/180);             % set arbritrary tumble of spacecraft [rad/s]
POLAR.L = [0;0;0];                                % sum of torques
POLAR.bias = randn(3,1)*1e-3;                     % gyro bias

POLAR.greg = datevec(datetime);
POLAR.JD = greg2jd(POLAR.greg);

% Sensor Specs
POLAR.sig_b = 1/300;            % magnetometer noise [Tesla]
POLAR.sig_s = 0.2*pi/180;         % sun sensor noise [radians]
POLAR.ss_boresight = [-1 0 0];    % sun-sensor boresight
POLAR.sig_u = 8.7e-7;             % gyro drift (in-run bias stability)
POLAR.sig_v = 3.8e-6;             % gyro noise

% Reaction Wheel Properties
POLAR.wheelaxis = eye(3);         % Alignment of each wheel, could support pyramid/NASA standard etc.
POLAR.wheelspeeds = [0;0;0];      % wheel speeds
POLAR.maxw_wheel = 850*[1;1;1];   % max wheel speed [rad/s]
POLAR.wheeltorques = [0;0;0];     % placeholder for dynamics
POLAR.maxT_wheel = 0.005*[1;1;1]; % max wheel torque [N m]
POLAR.wheel_J = 9.1e-6*[1;1;1];    % wheel inertia

% Magnetorquer Properties
POLAR.mag_ECI = [0;0;0];
POLAR.mag_moment = 0.2*[1;1;1]; % magnetic torquer max moment [A m^2]
POLAR.mag_axis = eye(3); 
POLAR.mag_duty = [0;0;0]; % duty cycle, between -1 and 1


%%% TARGET %%%
a = 6783;  %  altitude [km]
e = 1e-4;  %  eccentricity [unitless]
i = 51;    %  inclination  [degrees]
RAAN = 1.95;  %  RAAN         [degrees]
argp = -.6;  %  Arg. of perigee [degrees]
M = .1;    %  Mean Anomaly [radians]
TARGET.OE = [a,e,i,argp,RAAN,M];   % orbital elements of TARGET 
[r,v] = kep2rv(a,e,i,argp,RAAN,M);                                               
TARGET.r = r; TARGET.v = v;

TARGET.q = [1;0;0;0];
TARGET.omega = [.1;.2;.3]*pi/180;


% %%% Ground Station %%%
% a = ;
% e = ;
% i = ;
% RAAN = ;
% argp = ;
% M = ;
% GS.OE = [[a,e,i,argp,RAAN,M];
% [r,v] = kep2rv(a,e,i,RAAN,argp,M);
% GS.r = r; GS.v = v;

end

