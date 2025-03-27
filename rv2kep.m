function [OE] = rv2kep(r,v)
% This function calculates the orbital elements for a given a set of state
% vectors.
%
% This algorithm comes directly from "Fundamentals of Astrodynamics and
% Applications" by Vallado.

% Define standard gravitaitonal parameter:
mu = 398600.64;

% Calculate Angular Momentum:
H = cross(r,v);

% Calculate node vector:
n = cross([0 0 1],H);

% Calculate eccentricity:
evec = sqr((norm(v)^2-mu/norm(r))*r - dot(r,v)*v)/mu;
e    = norm(evec);
disp(evec)
% Calculate specific energy:
energy = (norm(v)^2)/2 - mu/norm(r);

% Calculate semi-major axis:
a = -mu/(2*energy);

% Calculate inclination:
inc = acosd(H(3)/norm(H));

% Calculate RAAN:
RAAN = acosd(n(1)/norm(n));
if n(2) < 0
   RAAN = 360-RAAN;
end

% Calculate Argument of Perigee:
argp = acos(dot(n,evec)/(norm(n)*e));
if evec(3) < 0
   argp = 360-argp;
end

% Calculate true anomaly:
nu = acosd(dot(evec,r)/(e*norm(r)));
if dot(r,v)<0
   nu = 360 - nu;
end

% Calculate true anomaly (Approximation):
M0 = nu - 2*e*sind(nu) + ((3/4)*(e^2) + (1/8)*e^4)*sind(2*nu) -...
     (1/3)*(e^3)*sind(3*nu) + (5/32)*(e^4)*sind(4*nu);

% Format output:
OE(1) = a;
OE(2) = e;
OE(3) = deg2rad(inc);
OE(4) = deg2rad(argp);
OE(5) = deg2rad(RAAN);
OE(6) = deg2rad(M0);

end