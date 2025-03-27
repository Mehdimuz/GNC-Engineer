function [Omega_bar] = Om_bar(w,dt)
%OM_BAR Summary of this function goes here
%   Detailed explanation goes here
 nw = norm(w); xi = sin(nw*dt/2)*w(:)/nw;
Omega_bar = eye(4)*cos(nw*dt/2) + [ -cpm(xi), xi; xi', 0];
end

