% lla = ecef2geodetic(r)
% By Benjamin Reifler 2017
% Modified by Chris Gnam 2019
% Adapted From Fundamentals of Spacecraft Attitude Determination and
%   Control by J. L. Crassidis and F. L. Markley
% Inputs:   r (km ECEF)
% Outputs:  lla = [lam phi h]' (deg & km geodetic)

function [lat, lon, alt] = ecef2wgs84(r) %#codegen

% calculate
a = 6378.1370;              % REquator
b = 6356.7523142;           % RPole
e = 0.081819190928905;      % sqrt(1 - b^2/a^2)
eps = 0.082094438036852;    % sqrt(a^2/b^2 - 1)
rho = hypot(r(1),r(2));
p = abs(r(3))/eps^2;
s = rho^2/(e^2*eps^2);
q = p^2 - b^2 + s;
u = p/sqrt(q);
v = b^2*u^2/q;
P = 27*v*s/q;
Q = (sqrt(P + 1) + sqrt(P))^(2/3);
t = (1 + Q + 1/Q)/6;
c = sqrt(u^2 - 1 + 2*t);
w = (c - u)/2;

% NEED TO INCLUDE abs() CALL HERE TO KEEP FROM GOING IMAGINARY WHEN r(3)~0
X = abs(sqrt(t^2 + v) - u*w - t/2 - 0.25);
d = sign(r(3))*sqrt(q)*(w + sqrt(X));
N = a*sqrt(1 + eps^2*d^2/b^2);
lat = asind((eps^2 + 1)*d/N);
lon = atan2(r(2), r(1))*180/pi;
alt = hypot(rho, r(3)) - a^2/N;

end
