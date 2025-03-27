function [FC, POLAR] = TrackingController(FC, POLAR, ref, meas)

quat_curr = FC.quat; w_curr = FC.w; wheel_speeds = POLAR.wheelspeeds;
quat_des = FC.desq; w_des = FC.desw; wheel_bias = FC.wheel_bias;
wheel_norms = POLAR.wheelaxis; wheel_inertia = POLAR.wheel_J;
mag_norms = POLAR.mag_axis; mag_field = meas.mag; 
integral_err = FC.integral_error; dt = FC.dt; 
kp = FC.kp; kd = FC.kd; ki = FC.ki; J = FC.J; L_prev = FC.L_prev;
%TRACKINGCONTROLLER PID framework to have camera track targets
%   Function requires quat_curr and w_curr to come from the MEKF and
%   quat_des and w_des to come from the MissionTargeter. Wheel speeds are
%   either tracked in real time, and wheel bias desired state for the
%   wheels (often 0,0,0). kp is proportional gain (torque to contribute
%   from attitude error), kw is derivative gain (torque from ang rate
%   error) and ki integral gain (torque to command due to accumulated error
%   over time)

cpm = @(r) [0 -r(3) r(2); r(3) 0 -r(1); -r(2) r(1) 0]; % cross product matrix
Xi  = @(q) [q(4)*eye(3)+cpm(q(1:3)); -q(1) -q(2) -q(3)];


%{ 
These methods didn't work but I saved it anyways

%eq parameters (EQ. FROM ...)              see blue-book section 7.3;
bigOm = @(w) [-cpm(w) w ; -w' 0];
q_dot = @(q) (1/2)*Xi(q)*w_curr;
qdd   = -1/4*dot(w_des,w_des)*quat_des;
qddot = 0; %assume zero
L1 = kp*eye(3);
L2 = kd*eye(3);

dq13 = Xi(quat_des).'*quat_curr;
dq4  = quat_curr.'*quat_des;
dW = ( w_des - w_curr );
integral_err = max(min(integral_err + dq13*dt + dW*dt,10),-10);

% L_des = cpm(w_curr)*J*w_curr+...
%         2*J*inv(Xi(quat_des)'*Xi(quat_curr))*...
%         ( ( (1/4)*(w_curr'*w_curr)*Xi(quat_des)') - ...
%           (Xi( 1/2*Xi(quat_des)*w_des)'*bigOm(w_curr) ) -...
%           (Xi(qdd).'- ...
%           (L1)*Xi(quat_des)'-...
%           L2 * ( (1/2)*Xi(quat_des)'*bigOm(w_curr)+...
%           Xi( 1/2*Xi(quat_des)*w_des ).') ) )* quat_des;
% L_des = L_des/100;

%}


%%%% PD portion
% Effectively use state feedback to convert the freely-rotating satellite
% into a spring-mass-damper analogue using wheel torque

quat_curr = Om_bar(w_curr, .5)*quat_curr; %propagating ahead slightly settles quicker
del_q13 = Xi(quat_des)'*quat_curr;
del_q4 = dot(quat_curr, quat_des);
del_w = w_curr - w_des;

L_bar = -kp*sign(del_q4)*del_q13 - kd*(1-del_q13.'*del_q13)*del_w;
FC.L_des = L_bar;

% This function sets wheel torques
POLAR = Wheel_Saturation(POLAR, L_bar, dt);

%%%% Dumping portion
% Uses the magnetorques to torque against the undesired momentum currently
% being stored in the wheels

delta_h = wheel_norms*diag(wheel_inertia)*(wheel_speeds - wheel_bias);
mag_mom = cross(delta_h, mag_field/norm(mag_field));
POLAR.mag_duty = min(max(mag_mom./POLAR.mag_moment,-1),1);


% This one is for the MEKF
FC.L_prev = POLAR.L;

end
