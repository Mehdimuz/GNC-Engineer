% r = geodetic2ecef(lla)
% By Benjamin Reifler
% Adapted From Fundamentals of Spacecraft Attitude Determination and
%   Control by J. L. Crassidis and F. L. Markley
% Inputs:   lla = [lam phi h]' (deg & km geodetic)
% Outputs:  r (km ECEF)

function r = wgs842ecef(lla) %# codegen

a = 6378.1370;          % REquator
%b = 6356.7523142;       % RPole
%f = 0.0033528129;       % (a - b)/a
e = 0.081819218070436;  % sqrt(f*(2 - f))
N = a/sqrt(1 - e^2*sind(lla(1))^2);
r = [(N + lla(3))*cosd(lla(1))*cosd(lla(2));...
    (N + lla(3))*cosd(lla(1))*sind(lla(2));...
    (N*(1 - e^2) + lla(3))*sind(lla(1))];
