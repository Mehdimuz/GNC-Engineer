function [POLAR, PWR] = PWR_Model(PWR, POLAR, FC, ref)
%PWR_MODEL tracks the solar panels, battery charge, voltage across each
%line, etc.
%   Detailed explanation goes here

A = q2a(POLAR.q);

% determine if TX radio should be on
jd = POLAR.JD - 2451545;
jdf = mod(jd, 1);
theta = mod(jdf + 0.7790572732640 + 0.00273781191135448*jd,1)*2*pi;
ECI2ECEF = [cos(theta) sin(theta) 0 ; -sin(theta) cos(theta) 0 ; 0 0 1];
GS_lat = FC.GS_lat; GS_lon = FC.GS_lon;
R_GS = 6371*ECI2ECEF'*[cosd(GS_lat)*cosd(GS_lon); cosd(GS_lat)*sind(GS_lon); sind(GS_lat)];
POLAR.elev_angle = asind(dot(POLAR.r-R_GS, R_GS)/norm(POLAR.r-R_GS)/norm(R_GS));
PWR.TX_on = POLAR.elev_angle > 10;
if PWR.TX_on
%keyboard
end
switch FC.mode
    case 'Mission'
        PWR.polarimeter_duty = ref.inSun;
        PWR.FC_duty = 1;
        PWR.ADCS_cdh_duty = 1;
        PWR.ADCS_mag_duty = 1;
        PWR.ADCS_sun_duty = ref.inSun;
        PWR.Radio_TX_duty = 0;
        PWR.Radio_RX_duty = 1;
        PWR.Radio_BL_duty = 0;
        PWR.Radio_BA_duty = 0;
        PWR.Radio_S_duty  = 0;
        PWR.GPS_duty = 0.2;
    case 'Sunpointing'
        PWR.polarimeter_duty = 0;
        PWR.FC_duty = 1;
        PWR.ADCS_cdh_duty = 1;
        PWR.ADCS_mag_duty = 1;
        PWR.ADCS_sun_duty = ref.inSun;
        PWR.Radio_TX_duty = PWR.TX_on*0.5;
        PWR.Radio_RX_duty = 1;
        PWR.Radio_S_duty  = PWR.TX_on*0.2;
        PWR.GPS_duty = 0.2;
    case 'Downlink'
        PWR.polarimeter_duty = 0;
        PWR.FC_duty = 1;
        PWR.ADCS_cdh_duty = 1;
        PWR.ADCS_mag_duty = 1;
        PWR.ADCS_sun_duty = ref.inSun;
        PWR.Radio_TX_duty = PWR.TX_on*0.5;
        PWR.Radio_RX_duty = 1;
        PWR.Radio_S_duty  = PWR.TX_on*0.2;
        PWR.GPS_duty = 0.2;
    case 'Detumble'
        PWR.polarimeter_duty = 0;
        PWR.FC_duty = 1;
        PWR.ADCS_cdh_duty = 1;
        PWR.ADCS_mag_duty = 1;
        PWR.ADCS_sun_duty = ref.inSun;
        PWR.Radio_TX_duty = 0;
        PWR.Radio_RX_duty = 1;
        PWR.Radio_S_duty  = 0;
        PWR.GPS_duty = 0.2;
    case 'Safe'
        PWR.polarimeter_duty = 0;
        PWR.FC_duty = 1;
        PWR.ADCS_cdh_duty = 0;
        PWR.ADCS_mag_duty = 0;
        PWR.ADCS_sun_duty = 0;
        PWR.Radio_TX_duty = 0;
        PWR.Radio_RX_duty = 1;
        PWR.Radio_S_duty  = 0;
        PWR.GPS_duty = 0.2;
    otherwise
        PWR.polarimeter_duty = 0;
        PWR.FC_duty = 0;
        PWR.ADCS_cdh_duty = 0;
        PWR.ADCS_mag_duty = 0;
        PWR.ADCS_sun_duty = 0;
        PWR.Radio_TX_duty = 0;
        PWR.Radio_RX_duty = 0;
        PWR.Radio_S_duty  = 0;
        PWR.GPS_duty = 0;
end

% rough linear fit from the data sheet
PWR.ADCS_wheels_i = (.12*abs(POLAR.wheel_J.*POLAR.wheelspeeds) ...
       -120*POLAR.wheeltorques.*POLAR.wheel_J.*POLAR.wheelspeeds)*...
       1000/PWR.eps_v_lines(end);

pc = [0 0 0 0];

pc(PWR.polarimeter_v) = pc(PWR.polarimeter_v) + PWR.polarimeter_i*PWR.polarimeter_duty;
pc(PWR.FC_v) = pc(PWR.FC_v) + PWR.FC_i*PWR.FC_duty;
pc(PWR.ADCS_cdh_v) = pc(PWR.ADCS_cdh_v) + PWR.ADCS_cdh_i*PWR.ADCS_cdh_duty;
pc(PWR.ADCS_mag_v) = pc(PWR.ADCS_mag_v) + PWR.ADCS_mag_i*PWR.ADCS_mag_duty;
pc(PWR.ADCS_sun_v) = pc(PWR.ADCS_sun_v) + PWR.ADCS_sun_i*PWR.ADCS_sun_duty;
pc(PWR.Radio_RX_v) = pc(PWR.Radio_RX_v) + PWR.Radio_RX_i*PWR.Radio_RX_duty;
pc(PWR.Radio_TX_v) = pc(PWR.Radio_TX_v) + PWR.Radio_TX_i*PWR.Radio_TX_duty;
pc(PWR.Radio_S_v) = pc(PWR.Radio_S_v) + PWR.Radio_S_i*PWR.Radio_S_duty;
pc(PWR.GPS_v) = pc(PWR.GPS_v) + PWR.GPS_i*PWR.GPS_duty;
pc(PWR.ADCS_mtq_v) = pc(PWR.ADCS_mtq_v) + sum(PWR.ADCS_mtq_i*abs(POLAR.mag_duty));
pc(PWR.ADCS_wheels_v) = pc(PWR.ADCS_wheels_v) + sum(PWR.ADCS_wheels_i);

PWR.eps_i_lines = pc;
PWR.consumed = dot(pc,PWR.eps_v_lines./PWR.eps_eff_lines);

PWR.sun_angs = max(PWR.panel_norms*(A*ref.sun),0);
PWR.generated = ref.inSun*1300*...
  dot(PWR.panel_areas.*PWR.panel_eff, PWR.sun_angs);
PWR.sun_angs = PWR.sun_angs + 0.01*randn(size(PWR.panel_eff));

PWR.bat_charge = PWR.bat_charge +FC.dt/3600*(PWR.generated - PWR.consumed);
PWR.bat_charge = max(min(PWR.bat_charge, PWR.bat_charge_max),0);

end

