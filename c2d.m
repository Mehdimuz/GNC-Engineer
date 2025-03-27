% [phi, gamma] = c2d_codegen(a, b, t)
% Continuous to discrete system conversion, adapted from MATLAB's c2d()
% Inputs:   a
%           b
%           t
% Outputs:  phi
%           gamma

function [phi, gamma] = c2d_codegen(a, b, t)

[~,n] = size(a);
[~,nb] = size(b);
s = expm([[a b]*t; zeros(nb,n+nb)]);
phi = s(1:n,1:n);
gamma = s(1:n,n+1:n+nb);
