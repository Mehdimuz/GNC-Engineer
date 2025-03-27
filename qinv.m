% This function computes the inverse quaternion.
% Andrew Dianetti
% 4 December 2014

function invq = qinv(q)

invq=[-q(1:3);q(4)];