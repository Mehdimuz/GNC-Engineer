% This function computes the Xi matrix for a quaternion.
% Andrew Dianetti
% 5 December 2014

% Required functions:
%   cpm.m - Calculates cross product matrix

function xi = Xi(q)
xi=[q(4)*eye(3)+cpm(q(1:3)); -q(1) -q(2) -q(3)]; 
