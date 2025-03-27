function [R_eci,V_eci] = kep2rv(a,e,i,omega,Omega,theta)
    %Created 7/11/2014 by Chris Shelton
    %Given the orbital elements this function will find the positition and
    %velocity in ECI (m)
    mu = 398600;

    % Define Orbit Parameters
    h = sqrt(mu.*a.*(1-e.*e));
    r = ((h.*h)./mu).*(1./(1+(e.*cos(theta))));
    % Calculate R,V in PQW coordinates
    R_pqw(:,1) = r.*cos(theta);
    R_pqw(:,2) = r.*sin(theta);
    R_pqw(:,3) = 0;
    V_pqw(:,1) = (mu./h).*-sin(theta);
    V_pqw(:,2) = (mu./h).* (e+cos(theta));
    V_pqw(:,3) = 0;

    i = i*pi/180;
    % Define the transfomation matrix from PQW to ECI
    a1 = (-sin(Omega).*cos(i).*sin(omega) + cos(Omega).*cos(omega));
    a2 = (-sin(Omega).*cos(i).*cos(omega) - cos(Omega).*sin(omega));          
    % a3 = (sin(i)*sin(omega));
    a4 = (cos(Omega).*cos(i).*sin(omega) + sin(Omega).*cos(omega));
    a5 = (cos(Omega).*cos(i).*cos(omega) - sin(Omega).*sin(omega));
    % a6 = (sin(i)*-cos(omega));
    a7 = (sin(omega).*sin(i));
    a8 = (cos(omega).*sin(i));
    % a9 = (cos(i));

    % Transform the vectors from PQW into ECI
    R_eci = [(a1.*R_pqw(:,1))+(a2.*R_pqw(:,2));(a4.*R_pqw(:,1)+a5.*R_pqw(:,2));(a7.*R_pqw(:,1) + a8.*R_pqw(:,2))];
    V_eci = [(a1.*V_pqw(:,1))+(a2.*V_pqw(:,2));(a4.*V_pqw(:,1)+a5.*V_pqw(:,2));(a7.*V_pqw(:,1) + a8.*V_pqw(:,2))];
end

