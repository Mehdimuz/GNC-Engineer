function [euler] = q2e(quat)
    A = q2a(quat);
    euler = a2e(A);
end