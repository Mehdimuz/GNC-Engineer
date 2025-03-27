function [angleError_degrees, q_error] = quatError(q1,q2)
% This function returns the angle, in degrees, that exists between two
% quaternions. Uses Shuster convention.
% Written by Nick DiGregorio

% First perform quaternion multiplication to get the error quaternion.
q_error = quatmult(q1, qinv(q2));


% Normalize the quaternion! It will fix all of your problems in GNC.
q_error = q_error / norm(q_error);


% Calculate a rotational angle error from the error quaternions, in
% degrees. These are the most useful metrics of error for humans to
% read. Reference: Dr. Crassidis lecture notes, lecture 17 Attitude
% Representations, page 21. This formulation makes use of the small
% angle approximation.
angleError_radians = 2 * q_error(1:3);
angleError_degrees = angleError_radians * 180/pi;

end