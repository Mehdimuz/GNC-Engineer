function FC = MissionTargeterSimplified(FC, POLAR, TARGET, ref, meas)
%   This function is like 'missiontargeter', but it assumes that the
%   boresight is directly along the body x and the desired spin (along with
%   the solar panel) is directly along the body z, which allows us to make
%   some helpful simplifications in the algorithm.
%


r = POLAR.r; v = POLAR.v;
target_R = TARGET.r; target_V = TARGET.v;

% [target_R, target_V] = kep2rv(targetOE(1),targetOE(2),targetOE(3),targetOE(4),targetOE(5),targetOE(6)); % convert target OE to pos, vel in ECI
LOS = target_R(:) - r(:);           % difference in POLAR/target's pos. vectors = LOS vector 
LOS_dot = target_V(:) - v(:);       % change rate of LOS 

body_x_des = LOS/norm(LOS); % desired "x" face is directly along LOS (for pointing)
body_z_des = cross(body_x_des, LOS_dot); % sets "z" face to track along spin axis 
% 
% if norm(LOS_dot)/norm(LOS) < 0.001
%     prop_target_R = target_R; prop_target_V = target_V;
%     prop_polar_R  = r;        prop_polar_V  = v;
%     dxo = zeros(6,1);
%     dxp = zeros(6,1);
%     temp_ = 0;
%     while temp_ < 0.001
%         % propagate target's orbit forward by 1 second
%         k1o = sub_dynamics(prop_target_R,prop_target_V);
%         k2o = sub_dynamics(prop_target_R+k1o(1:3)/2,prop_target_V+k1o(4:6)/2);
%         k3o = sub_dynamics(prop_target_R+k2o(1:3)/2,prop_target_V+k2o(4:6)/2);
%         k4o = sub_dynamics(prop_target_R+k3o(1:3),prop_target_V+k3o(4:6));
%         dxo = dxo + k1o + 2*k2o + 2*k3o + k4o;
%         prop_target_R = target_R + dxo(1:3);
%         prop_target_V = target_V + dxo(4:6);
%         % now polar
%         k1p = sub_dynamics(prop_polar_R,prop_polar_V);
%         k2p = sub_dynamics(prop_polar_R+k1p(1:3)/2,prop_polar_V+k1p(4:6)/2);
%         k3p = sub_dynamics(prop_polar_R+k2p(1:3)/2,prop_polar_V+k2p(4:6)/2);
%         k4p = sub_dynamics(prop_polar_R+k3p(1:3),prop_polar_V+k3p(4:6));
%         dxp = dxp + k1p + 2*k2p + 2*k3p + k4p;
%         prop_polar_R = r + dxp(1:3);
%         prop_polar_V = v + dxp(4:6);
%         temp_ = norm(prop_target_V-prop_polar_V)/norm(prop_target_R-prop_polar_R);
%     end
%     body_z_des = cross(body_x_des, prop_target_V - prop_polar_V);
% end

% construct desired attitude matrix
body_z = LOS;
body_y = cross(body_z, LOS_dot);
body_x = cross(body_y, body_z);
A = [body_x body_y body_z];
A = (A./vecnorm(A))';

if isempty(FC.desq)
    FC.desq = a2q(A);
    FC.desw = [0;0;1].*norm(LOS_dot)/norm(LOS);
else % low pass-filter, to keep from desired attitude jittering
    a = 0.5; % weight of the dynamics
    desq = a*(Om_bar(FC.desw,FC.dt)*FC.desq) + (1-a)*a2q(A);
    FC.desq = desq/norm(desq)*sign(desq(4));
    FC.desw = a*FC.desw + (1-a)*[0;0;1]*norm(LOS_dot)/norm(LOS);
end
FC.los = LOS;
FC.los_dot = LOS_dot;

end

function dx = sub_dynamics(r,v)
mu = 398400;

% derivative of each state
dr = v;
dv = -mu/norm(r)^3*r;

dx = [dr; dv];
end

