% By Benjamin Reifler
% Note that this is for direction vectors, to use position vectors,
%   subtract the reference point first
% Inputs:   ned (3x1 vector in NED)
%           lat (deg)
%           lon (deg)
%           JD (Julian date)
%           coefs (1600x17,1275x17,66x11)
% Outputs:  A (DCM)

function [v_eci] = ned2eci(v_ned, lat, lon, eci2ecef) %# codegen
    % NED -> ECEF rotation matrix
    cosl = cosd(lon);
    sinl = sind(lon);
    cosp = cosd(lat);
    sinp = sind(lat);
    A = [-sinp*cosl -sinl  -cosp*cosl;
         -sinp*sinl  cosl  -cosp*sinl;
            cosp      0      -sinp];

    % NED -> ECEF -> ECI rotation
    v_eci = eci2ecef'*A*v_ned;
end
