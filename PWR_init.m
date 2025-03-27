function PWR = PWR_init(POLAR)

% These are rough estimates
PWR.panel_norms = [0 0 1; 0 1 0; 0 -1 0];
PWR.panel_areas = [0.1*0.3; 0.1*0.3; 0.1*0.3];
PWR.panel_eff   = [0.2; 0.2; 0.2];
PWR.max_gen = 8.4;

PWR.bat_charge = 20;
PWR.bat_charge_max = 30;
PWR.vat_voltage = 16;
PWR.eps_v_lines = [3.3, 5, 12, 8];
PWR.eps_i_lines = [0,0,0,0];
PWR.eps_eff_lines = [0.88 0.91 0.93 1];

% component power draw
% _v is which power line it's on
% _i is current draw in amps
PWR.polarimeter_v = 3;  PWR.polarimeter_i = 0.21;
PWR.FC_v = 2;           PWR.FC_i = 0.25; % Flight Computer
PWR.ADCS_cdh_v = 1;     PWR.ADCS_cdh_i = 0.07; % CubeComputer
PWR.ADCS_mag_v = 1;     PWR.ADCS_mag_i = 0.015;% CubeMag
PWR.ADCS_sun_v = 1;     PWR.ADCS_sun_i = 0.028;% CubeSense Sun
PWR.Radio_TX_v = 2;     PWR.Radio_TX_i = 0.95; % PULSAR UTRXC
PWR.Radio_RX_v = 2;     PWR.Radio_RX_i = 0.04; % UHF Radio
PWR.Radio_S_v  = 4;     PWR.Radio_S_i  = 1.5; % ISIS S-band TX
PWR.GPS_v = 1;          PWR.GPS_i = 0.01; % U-blox m8n

PWR.polarimeter_duty = 0;
PWR.FC_duty = 0;
PWR.ADCS_cdh_duty = 0;
PWR.ADCS_mag_duty = 0;
PWR.ADCS_sun_duty = 0;
PWR.Radio_TX_duty = 0;
PWR.Radio_RX_duty = 0;
PWR.Radio_BL_duty = 0;
PWR.Radio_BA_duty = 0;
PWR.Radio_S_duty  = 0;
PWR.GPS_duty = 0;


% ADCS actuators have highly variable duty cycle
PWR.ADCS_wheels_v = 4; % only thing on Vbat line
PWR.ADCS_mtq_v = 2; PWR.ADCS_mtq_i = 0.1;

% delta at each timestep
PWR.generated = 0;
PWR.consumed = 0;

max_gen = 0;
for theta = 1:360
    for phi = -90:90
        r = [cosd(theta)*cosd(phi); sind(theta)*cosd(phi); sind(phi)];
        r = r/norm(r);
        gen =dot(PWR.panel_areas.*PWR.panel_eff, max(PWR.panel_norms*r,0));
        if gen > max_gen    
            max_gen = gen;
            PWR.des_sun = r;
        end
    end
end

end

