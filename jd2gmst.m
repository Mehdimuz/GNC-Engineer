function [GMST] = JD2GMST(JD)
    T_UT1 = (JD - 2451545.0)/36525.0;

    GMST = - 6.2e-6 * T_UT1 * T_UT1 * T_UT1 + 0.093104 * T_UT1 * T_UT1  ...
           + (876600.0 * 3600.0 + 8640184.812866) * T_UT1 + 67310.54841;

    GMST = rem(deg2rad(GMST)/240.0, 2*pi);

    if (GMST < 0)
        GMST = GMST + 2*pi;
    end
end

