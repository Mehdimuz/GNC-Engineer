function [desq, desw, LOS, LOS_dot] = MissionTargeter(r, v, boresight, J_max, targetOE)
%TARGETER: determines necessary quaternion and angular rate for POLAR pointing
%   INPUTS -
%   q : current satellite quaternion (4x1)
%   w : current satellite angular rate in rad/s (3x1)
%   r : current satellite position in km in ECI (3x1)
%   v : current satellite velocity in km/s in ECI (3x1)
%   boresight: vector of payload boresight in body frame (3x1)
%   J_max : vector of max principal inertia (3x1)
%   targetOE : orbital elements of the targeter 
%             (this should come from the TLE)
%         [ semimajor-axis (km);
%                   ecentricity;
%              inclination(rad);
%          right-ascension(rad);
%      argument of perigee(rad);
%             mean anomaly(rad)];
%
%   OUTPUTS - 
%   desq : desired quaternion (4x1)
%   desw : desired angular rate (rad/s) (3x1)

[target_R, target_V] = kep2rv(targetOE(1),targetOE(2),targetOE(3),targetOE(4),targetOE(5),targetOE(6));
LOS = target_R - r;
LOS_dot = target_V - v;

% first, perform a rotation to align the boresight vector with the line of
% sight vector
boresight = boresight/norm(boresight);
e1 = cross(boresight,LOS);
e1 = e1/norm(e1);
theta1 = acos(dot(boresight, LOS)/(norm(boresight)*norm(LOS)));
q1 = [e1*sin(theta1/2); cos(theta1/2)];

% This frame only has a one axis constraint, and is free to rotate about
% the payload boresight. The final axes are determined so that the
% desired angular rate falls about the desired spin axis.
A = q2a(q1);

body_LOS_dot = A*LOS_dot/norm(LOS);
% solving the problem where there is nearly no angular rate in the target
% tracking causing instability, and propagating forward
if norm(body_LOS_dot) < 0.001
    prop_target_R = target_R; prop_target_V = target_V;
    dx = zeros(6,1);
    i = 0;
    while norm(body_LOS_dot) < 0.001
        % propagate targets orbit forward by 1 second
        k1 = sub_dynamics(prop_target_R,prop_target_V);
        k2 = sub_dynamics(prop_target_R+k1(1:3)/2,prop_target_V+k1(4:6)/2);
        k3 = sub_dynamics(prop_target_R+k2(1:3)/2,prop_target_V+k2(4:6)/2);
        k4 = sub_dynamics(prop_target_R+k3(1:3),prop_target_V+k3(4:6));
        dx = dx + k1 + 2*k2 + 2*k3 + k4;
        prop_target_R = prop_target_R + dx(1:3);
        prop_target_V = prop_target_V + dx(4:6);
        body_LOS_dot = A*(LOS_dot + dx(1:3))/norm(LOS);
        i = i + 1;
    end
end

e2 = cross(boresight, body_LOS_dot);
e2 = e2/norm(e2);

% make the spin axis perpendicular to the boresight
s = cross(boresight, cross(J_max, boresight));
s = s/norm(s);

theta2 = acos(dot(s, e2));
q2 = [boresight*sin(theta2/2); cos(theta2/2)];

cpm = @(r) [0 -r(3) r(2); r(3) 0 -r(1); -r(2) r(1) 0];
Psi = @(q) [-cpm(q)+eye(3)*q(4);-q(1) -q(2) -q(3)];
q_cross=@(q)[Psi(q) q];
desq = q_cross(q2) *q1;
p = LOS_dot/norm(LOS);
desw = s*norm(p);

end

function dx = sub_dynamics(r,v)
mu = 398400;

% derivative of each state
dr = v;
dv = -mu/norm(r)^3*r;

dx = [dr; dv];
end