function [rho] = AtmosDensModel(alt)

% Atmospheric Density Model Based on Vallado's Code

if alt >= 800
   H = 130.8;
   rho_i = 4.262*10^-17;
   rho = rho_i*exp((800-alt)/H);

elseif alt >= 700
   H = 105.3; 
   rho_i = 1.216*10^-16; 
   rho = rho_i*exp((700-alt)/H);

elseif alt >= 600
    H = 91;
    rho_i = 3.818*10^-16;
    rho = rho_i*exp((600-alt)/H);

elseif alt >= 500
    H = 81.9;
    rho_i = 1.316*10^-15;
    rho = rho_i*exp((500-alt)/H);

elseif alt >= 400
    H = 73.2;
    rho_i = 5.192*10^-15;
    rho = rho_i*exp((400-alt)/H);

elseif alt >= 300 
    H = 61.2;
    rho_i = 2.653*10^-14;
    rho = rho_i*exp((300-alt)/H);

elseif alt >= 250
    H = 52.6;
    rho_i = 7.316*10^-14;
    rho = rho_i*exp((250-alt)/H);

elseif alt >= 200
    H = 40.8;
    rho_i = 2.706*10^-13;
    rho = rho_i*exp((200-alt)/H);

elseif alt >= 150
    H = 24.1;
    rho_i = 2.141*10^-12;
    rho = rho_i*exp((150-alt)/H);

elseif alt >= 130
    H = 16.1;
    rho_i = 8.484*10^-12;
    rho = rho_i*exp((130-alt)/H);

else
    H = 8.06;
    rho_i = 9.661*10^-11;
    rho = rho_i*exp((110-alt)/H);

    
end


end