function [reference, measurements,POLAR] = read_sensors(POLAR, PWR, dt, coefs)

A = q2a(POLAR.q); % ECI2BODY

jd = POLAR.JD - 2451545;
jdf = mod(jd, 1);

% Earth rotation angle
% this doesn't account for Earth's nutation, but that's a lot more
% precision than we need
theta = mod(jdf + 0.7790572732640 + 0.00273781191135448*jd,1)*2*pi;
ECI2ECEF = [cos(theta) sin(theta) 0 ; -sin(theta) cos(theta) 0 ; 0 0 1];


% Generate magnetometer readings
r_ecef = ECI2ECEF*POLAR.r;
[lat,lon,alt] = ecef2wgs84(r_ecef);
reference.mag = igrf(coefs, lat, lon, alt).';
measurements.mag = A*reference.mag + 1000*POLAR.sig_b*randn(3,1);  % magnetometer measurement                

% Generate sun-sensor reading

spring_jd = POLAR.JD - greg2jd([POLAR.greg(1), 3, 21, 0, 0, 0]);
spring_jdf = mod(spring_jd/365, 1);
reference.sun = [1;0;0]*cos(2*pi*spring_jdf) + [0;cosd(23.5);sind(23.5)]*sin(2*pi*spring_jdf);

% check if in eclipse
reference.inSun = true;
if dot(POLAR.r, reference.sun) > 0
    proj = (reference.sun.'*POLAR.r)*reference.sun;
    if norm(POLAR.r-proj) < 6371 
        reference.inSun = false;
    end
end

measurements.sun = nan(3,1);
measurements.sunAvail = false;
if reference.inSun
    if dot(POLAR.ss_boresight, A*reference.sun) > 0
        measurements.sun = A*reference.sun + POLAR.sig_s*randn(3,1);
        measurements.sunAvail = true;
    end
    S_angs = max(PWR.panel_norms*(A*reference.sun),0);
    reference.panel_norms = PWR.panel_norms(S_angs>.1,:);
    measurements.panel_angs = max(reference.panel_norms*(A*reference.sun),0);
    measurements.panel_angs = measurements.panel_angs + 0.01*randn(1,sum(S_angs>.1));
end

% gyroscope readings
measurements.omega = POLAR.omega + POLAR.sig_v*randn(3,1) + POLAR.bias;
POLAR.lat = lat;
POLAR.lon = lon;
end