function FC = MEKF(FC, ref, meas)
%MEKF THE MULTIPLICATIVE EXTENDED KALMAN FILTER

% extract from the structs
quat = FC.quat; omega = meas.omega; bias = FC.bias;
P = FC.P; Q = FC.Q; dt = FC.dt;
meas_mag = meas.mag; ref_mag = ref.mag; sig_mag = FC.sig_b;
meas_sun = meas.sun; ref_sun = ref.sun; sig_sun = FC.sig_s;
J = FC.J; L = FC.L_prev; omega_prev = FC.w;

%   quat - 4x1 a priori estimated quaternion
%   omega- 3x1 angular rate (rad/s)
%   bias - 3x1 angular rate bias (rad/s)
%   P    - state covariance (6x6)
%   Q    - process noise covariance (6x6)
%   dt   - timestep (1x1) (seconds)
%   
%   meas_- vectors coming from the sensors(3xn)
%   ref_ - vectors coming from environment models(3xn)
%   sig_ - uncertainty with sensors measurement(1xn or 3xn)
%
%   sig_u- gyro rate uncertainty (rad/s)
%   sig_v- gyro bias random walk (idr)
%   ^ The above two are covered in the Q matrix!
%   
%   J    - moment of inertia (kg m^2) -- currently unused
%   L    - sum of control torques (N m) -- currently unused
% 


% first, propagate the attitude to the current timestep
dqdt = @(q,w) 1/2*Xi(q)*w;
dwdt = @(w) J\(cross(w,J*w) + L);
W = omega_prev;

% here, it's propagating the attitude by interpolating the current and
% previous angular rate estimate. the commented out lines on the right and
% to propagate the angular rate through its equatons of motion, and is the
% only portion you need the moment of inertia and torques in attitude
% estimation
kq1 = dqdt(quat,W); % kw1 = dwdt(W);
kq2 = dqdt(quat+kq1*dt/2,(W+omega-bias)/2);% kw2 = dwdt((W+omega-bias)/2);
kq3 = dqdt(quat+kq2*dt/2,(W+omega-bias)/2); % kw3 = dwdt((W+omega-bias)/2);
kq4 = dqdt(quat+kq3*dt, omega-bias); % kw4 = dwdt(W+kw3*dt);
quat = quat + dt/6*(kq1 + 2*kq2 + 2*kq3 + kq4);
% nw = norm(omega-bias); psi = sin(nw*dt/2)*(omega-bias)/nw;
% quat = [eye(4)*cos(nw) + [-cpm(psi) psi; -psi.' 0] ]*quat;
quat = quat/norm(quat); % always renormalize your quaternion

% the state is a 3x1 of your attitude rotation error, and and 3x1 of your
% gyro bias. F is the linearized jacobian of the state (d/dt x = F*x)
F = [-cpm(omega-bias) eye(3); zeros(3,6)];
% Phi is the state transition matrix, estimated through a series approach.
% There is an analytical solution, this just saves some text & time
Phi = eye(6) + F*dt + F^2*dt^2/2 + F^3*dt^3/6;
% G is a matrix which expresses how noise in the process affects the state
G = diag([-1 -1 -1 1 1 1]);
% Propagate the covariance
P =  Phi*P*Phi'+ G*Q*G.';

% Create the attitude matrix and initialize the state
A = q2a(quat);
x = zeros(6,1);

% For each valid sensor measurement, the process is followed:
% normalize reference and measured vectors
% construct measurement covariance matrix R
% construct observation matrix H
% compute Kalman Gain matrix K
% compute a posteriori state through innovations process ( K*(y-H*x) )
% compute a posteriori state covariance P

if all(~isnan(meas_mag))
    meas_mag = meas_mag/norm(meas_mag);
    ref_mag = ref_mag/norm(ref_mag);
    R = diag([1 1 1].*sig_mag.^2);
    H = [cpm(A*ref_mag) zeros(3)];
    K = P*H.'/(H*P*H.'+R);
    y = [meas_mag];
    x = x + K*(y-H*x);
    P = P - K*H*P;
end

if meas.sunAvail
    % Fine Sun Sensor
    meas_sun = meas_sun/norm(meas_sun);
    ref_sun = ref_sun/norm(ref_sun);
    R = diag([1 1 1].*sig_sun.^2);
    H = [cpm(A*ref_sun) zeros(3)];
    K = P*H.'/(H*P*H.'+R);
    y = [meas_sun];
    x = x + K*(y-H*x);
    P = P - K*H*P;

    % Coarse Sun Sensors on Panels
    for jj = 1:length(meas.panel_angs)
        y = meas.panel_angs(jj);
        H = [ref.panel_norms(jj,:)*cpm(ref.sun), zeros(1,3)];
        K = P*H.'/(H*P*H.' + 0.1);
        x = x + K*(y-H*x);
        P = P - K*H*P;
    end
end

% Startracker is omitted, if you want it do this
%{ 
if all(~isnan(meas_stars))
    for i = 1:num_stars
    meas_stars(:,i) = meas_stars(:,i)/norm(meas_stars(:,i));
    ref_stars(:,i) = ref_stars(:,i)/norm(ref_stars(:,i));
    R = diag([1 1 1].*sig_star.^2);
    H = [cpm(A*ref_stars(:,i)) zeros(3)];
    K = P*H.'/(H*P*H.'+R);
    y = [meas_stars(:,i)];
    x = x + K*(y-H*x);
    P = P - K*H*P;
    end
end
%}

da = x(1:3); % extract error vector
quat = quat + 1/2*Xi(quat)*da; % update quaternion
FC.quat = quat/norm(quat)*sign(quat(4));
FC.bias = bias - x(4:6); % update bias
FC.w = omega - bias;
FC.P = P;

end
