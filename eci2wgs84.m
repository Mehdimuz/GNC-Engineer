% lla = eci2geodetic(r,JD,coefs)
% By Benjamin Reifler
% Inputs:   r (km ECI)
%           JD (Julian date)
%           coefs (1600x17,1275x17,66x11)
% Outputs:  lla = [lat lon h]' (deg & km geodetic)

function lla = eci2wgs84(r, eci2ecef)
    p = eci2ecef*r;

    % convert ECEF to geodetic
    lla = ecef2wgs84(p);
end