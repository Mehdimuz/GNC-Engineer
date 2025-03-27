function [JD] = greg2jd(year,month,day,hour,minute,sec)

if nargin == 1
    sec = year(6);
    minute = year(5);
    hour = year(4);
    day = year(3);
    month = year(2);
    year = year(1);
end
    % Now break the conversion up into different terms 
    term2 = floor( 7*(year + floor( (month+9) /12) ) /4 );
    term3 = floor( 275*month / 9);

    % Perform the calculation
    JD = 367*year - term2 + term3 + day + 1721013.5 +...
         hour/24 + minute/1440 + sec/86400;
end