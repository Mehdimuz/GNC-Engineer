% r = geodetic2eci(lla,JD,xCoefs,yCoefs,sCoefs)
% By Benjamin Reifler
% Inputs:   lla = [lat lon h]' (deg & km geodetic)
%           JD (Julian date)
%           coefs (1600x17,1275x17,66x11)
% Outputs:  r (km ECI)

function r = wgs842eci(lla, eci2ecef) %# codegen

% convert geodetic to ECEF
p = geodetic2ecef(lla);

r = eci2ecef'*p;
