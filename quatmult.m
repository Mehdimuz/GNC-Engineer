% This function returns the product of two quaternions.  Uses Shuster
% convention instead of traditional Hamiltonian convention.

% Required functions:
%     cpm.m - Computes cross product matrix

% Andrew Dianetti
% 2 December 2014

function prod = quatmult(q,r)

phi=[q(4)*eye(3)-cpm(q(1:3)); -q(1:3)'];
prod=[phi q]*r;
prod=prod/norm(prod); %Normalize the quaternion to eliminate machine precision issues
