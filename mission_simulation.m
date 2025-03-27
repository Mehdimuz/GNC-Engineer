clearvars -except MC LOGLOG; close all; clc;


[POLAR, TARGET] = mission_init;     % Initiate POLAR parameters
[FC, coefs, stars, mags] = GNC_init(POLAR);
PWR = PWR_init;

endtime = 10*60*60;   % PDR EDIT --> should run for ~8 hours (crassidis)
tspan = 0:FC.dt:endtime;            % Timespan of sim [s]


if exist('MC', 'var')
    monte_carlize
    fprintf("Run number: %d\n", MC.runs)
end

% Initialization for logging over time
Lsim = length(tspan);               % length of sim

% Dynamics-based fields
LOG.position = zeros(3,Lsim);       % position (true)
LOG.velocity = zeros(3,Lsim);		% velocity (true)
LOG.trueq = zeros(4,Lsim);          % true quaternion
LOG.omega = zeros(3,Lsim);          % ang. rate (true)
LOG.H = zeros(3,Lsim);              % body ang. momentum (true)
LOG.h = zeros(3,Lsim); 			    % wheels ang. momentum (N m s)
LOG.Lw = zeros(3,Lsim);             % wheel torques (N m)

% Sensor Logs
LOG.inSun = zeros(1,Lsim);		    % true in sun, false in eclipse
LOG.sunAvail = zeros(1,Lsim);		% true if sunsensor sees sun
LOG.bref = zeros(3,Lsim);           % mag vector in ECI
LOG.sref = zeros(3,Lsim);           % sun vector in ECI
LOG.bbod = zeros(3,Lsim);           % mag vector in body
LOG.sbod = zeros(3,Lsim);           % sun vector in body
LOG.W = zeros(3,Lsim);              % measured angular rate in body
LOG.bias = zeros(3,Lsim);		    % true gyro bias

% Estimation-based fields
LOG.estq =  zeros(4,Lsim);          % estimated quaternion
LOG.w = zeros(3,Lsim);              % estimated angular rate
LOG.wb = zeros(3,Lsim);			    % gyro bias
LOG.sig3 = zeros(6,Lsim); 		    % 3-sigma bounds

% Controller-based fields
LOG.q_des = zeros(4,Lsim);		    % Desired quaternion for pointing
LOG.w_des = zeros(3,Lsim);		    % Desired angular rate for pointing
LOG.los = zeros(3,Lsim); 		    % Line-of-sight (POLAR -> Target)
LOG.los_d = zeros(3,Lsim);		    % d/dt(Line-of-sight)
LOG.L_des = zeros(3,Lsim);
LOG.L = zeros(3,Lsim); 			    % commanded wheel torques
Log.mu = zeros(3,Lsim);			    % commanded magnetic moments

% Error-based fields
LOG.qerr =  zeros(4,Lsim);          % error quaternion (trueq->q_est)
LOG.qerr_est = zeros(4,Lsim);		% error quaternion (q_des->q_est)
LOG.qerr_con = zeros(4,Lsim); 		% error quaternion (q_des->trueq)

% Power-based fields
LOG.pwr_gen = zeros(1,Lsim);        % Power generation by panels (Watts)
LOG.pwr_con = zeros(1,Lsim);        % Power consumed by satellite (Watts)
LOG.bat = zeros(1,Lsim);            % Charge of Batter (Watt hours)
LOG.wheel_i = zeros(3,Lsim);        % Current from Reaction Wheels (Amps)

% Pre Allocating LAT & Lon
LOG.lat = zeros(1, Lsim);           % Lat of GS
LOG.lon = zeros(1, Lsim);           % Lon of GS


%%% BEGIN SIMULATION %%%

tic
for i = 1:length(tspan)

    %%% Propagate Dynamics %%%
    [POLAR, TARGET] = dynamics(FC.dt, POLAR, TARGET);
    POLAR.JD = greg2jd(POLAR.greg +[0 0 0 0 0 i*FC.dt]); % this is done for better precision
    A = q2a(POLAR.q);
 
    %%% Read Sensors %%%
    [ref,meas,POLAR] = read_sensors(POLAR, PWR, FC.dt, coefs);
    POLAR.mag_ECI = ref.mag*1e-9;


    % POWER Model
    [POLAR, PWR] = PWR_Model(PWR, POLAR, FC, ref);

    %%% Run Algos %%%
    % attitude estimation
    FC = MEKF(FC, ref, meas);

    % Controls

    switch FC.mode
        case "Mission"
            FC = MissionTargeterSimplified(FC, POLAR, TARGET, ref, meas);
            [FC, POLAR] = TrackingController(FC, POLAR, ref, meas);
        case "Sunpointing"
            [FC, POLAR] = SunPointingController(FC, POLAR, PWR, ref, meas);
        case "Downlink"
            FC = DownlinkTargeter(FC, POLAR, ref, meas);
            [FC, POLAR] = TrackingController(FC, POLAR, ref, meas);
        case "Safe"

    end
    % POLAR.L = FC.L_des;
    POLAR.L = POLAR.wheeltorques + cross(POLAR.mag_moment.*POLAR.mag_duty,A*POLAR.mag_ECI);

    %%% Log stuff %%%
    % see above for LOG field definitions
    LOG.position(:,i) = POLAR.r;
    LOG.velocity(:,i) = POLAR.v;
    LOG.trueq(:,i) = POLAR.q;
    LOG.omega(:,i) = POLAR.omega;
    LOG.H(:,i) = POLAR.I*POLAR.omega;
    LOG.h(:,i) = POLAR.wheelaxis*diag(POLAR.wheel_J)*POLAR.wheelspeeds;
    LOG.Lw(:,i) = POLAR.wheeltorques;


    LOG.estq(:,i) = FC.quat;
    LOG.w(:,i) = FC.w;
    LOG.estbias(:,i) = FC.bias;
    LOG.bias(:,i) = POLAR.bias;
    LOG.sig3(:,i) = 3*sqrt(diag(FC.P));


    LOG.inSun(i) = ref.inSun;
    LOG.sunAvail(i) = meas.sunAvail;
    LOG.bref(:,i) = ref.mag;
    LOG.sref(:,i) = ref.sun;
    LOG.bbod(:,i) = meas.mag;
    LOG.sbod(:,i) = meas.sun;
    LOG.W(:,i) = meas.omega;

    LOG.q_des(:,i) = FC.desq;
    LOG.w_des(:,i) = FC.desw;
    LOG.los(:,i) = FC.los;
    LOG.los_d(:,i) = FC.los_dot;
    LOG.L(:,i) = POLAR.L;
    LOG.L_des(:,i) = FC.L_des;
    LOG.mu(:,i) = POLAR.mag_duty.*POLAR.mag_moment;

    LOG.qerr(:,i) = q_cross(qinv(POLAR.q))*FC.quat;
    LOG.q_err_est(:,i) = q_cross(qinv(FC.desq))*FC.quat;
    LOG.q_err_con(:,i) = q_cross(qinv(FC.desq))*POLAR.q;

    LOG.pwr_gen(i) = PWR.generated;
    LOG.pwr_con(i) = PWR.consumed;
    LOG.bat(i) = PWR.bat_charge;
    LOG.wheel_i(:,i) = PWR.ADCS_wheels_i;
    
    LOG.lat(i) = POLAR.lat;
    LOG.lon(i) = POLAR.lon;
    LOG.elev(i)= POLAR.elev_angle;
    LOG.TX(i) = PWR.TX_on;
end
toc

if ~exist('MC','var')
    doAnimation = false;
    Plotting_Results
end