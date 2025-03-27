function [unix] = JD2Unix(JD)
    unix = (JD - 2440587.5)*86400;
end