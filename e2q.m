function [quat] = e2q(Rx,Ry,Rz,rad)
    if nargin == 1
        Ry = Rx(2);
        Rz = Rx(3);
        Rx = Rx(1);
        rad = false;
    elseif nargin == 3
        rad = false;
    elseif nargin == 4
        % Do nothing
    else
        error('Incorrect number of input arguments')
    end
    A = e2a(Rx,Ry,Rz,rad);
    quat = a2q(A);
end