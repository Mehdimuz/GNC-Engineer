function [FC, coefs, stars, mags] = GNC_init(POLAR)
% Struct "FC" is supposed to mean "Flight Computer

FC.quat = POLAR.q+randn(4,1)*1e-3; FC.quat = FC.quat/norm(FC.quat);% a priori est. quaternion (todo: implement Davenport's Q method for this)
FC.w = POLAR.omega+POLAR.bias;                                     % angular rate of spacecraft including biases
FC.dt = 5;                                                         % timestep [s]
FC.bias = [0;0;0];                                                 % 

FC.J = POLAR.I;

% Sensor Specs
FC.sig_b = 10/300;    % magnetometer noise
FC.sig_s = 0.2*pi/180;              % sun sensor noise
FC.sig_st = NaN ;                   % star tracker noise 
FC.sig_u = 8.7e-6;                  % gyro drift (in-run bias stability)
FC.sig_v = 3.8e-5;                  % gyro noise (ARW)

% MEKF inits
FC.P = blkdiag(eye(3)*.01, 1e-6*eye(3));                     % state Covariance Matrix 
FC.Q = [(FC.sig_v^2*FC.dt + 1/3*FC.sig_u^2*FC.dt^3)*eye(3), 1/2*(FC.sig_u^2*FC.dt^2)*eye(3);     % process noise covariance
                    1/2*(FC.sig_u^2*FC.dt^2)*eye(3),            FC.sig_u^2*FC.dt*eye(3)];
FC.desq = [];
FC.desw = [];

% Controller gains/inits
FC.kp = 0.0005; % proportional gain
FC.kd = 0.0028; % derivative gain
FC.ki = 0; % integral gain
FC.integral_error = [0;0;0];
FC.L_prev = [0;0;0];
FC.wheel_bias = [0;0;0];

load("igrf_coefs.mat")   % load IGRF magnetic field model
load("brightest_stars.mat") % all stars brighter than magnitude 5

% Ground station
FC.GS_lat = 43;
FC.GS_lon = -78.3;

% Mode Select
% acceptable modes are Detumble, Checkout, Mission, Sunpointing, Safe
FC.mode = "Sunpointing";
