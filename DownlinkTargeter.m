function FC = DownlinkTargeter(FC, POLAR, GS, ref, meas)

r = POLAR.r; v = POLAR.v;

lat = FC.GS_lat; lon = FC.GS_lon;
e = sqrt(1-6356.752^2/6378.137^2);
R_mag = 6378.137*(1-e^2)/(1+e*cosd(lat));
R_ECEF = R_mag*[cosd(lat)*cosd(lon);
                cosd(lat)*sind(lon);
                          sind(lat)];

jd = POLAR.JD - 2451545;
jdf = mod(jd, 1);
theta = mod(jdf + 0.7790572732640 + 0.00273781191135448*jd,1)*2*pi;
ECI2ECEF = [cos(theta) sin(theta) 0 ; -sin(theta) cos(theta) 0 ; 0 0 1];

R_ECI = ECI2ECEF'*R_ECEF;
V_ECI = cross([0; 0; 2*pi/86344], R_ECI);

LOS = R_ECI-r;
LOS_dot = V_ECI-v;


% r = POLAR.r; v = POLAR.v;
% target_R = TARGET.r; target_V = TARGET.v;
% 
% % [target_R, target_V] = kep2rv(targetOE(1),targetOE(2),targetOE(3),targetOE(4),targetOE(5),targetOE(6)); % convert target OE to pos, vel in ECI
% LOS = target_R(:) - r(:);           % difference in POLAR/target's pos. vectors = LOS vector 
% LOS_dot = target_V(:) - v(:);       % change rate of LOS 
% 
% body_x_des = LOS/norm(LOS); % desired "x" face is directly along LOS (for pointing)
% body_z_des = cross(body_x_des, LOS_dot); % sets "z" face to track along spin axis 

body_y = LOS;
body_z = cross(body_y, LOS_dot);
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