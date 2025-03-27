function [A] = e2a(Rx,Ry,Rz,rad)
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
    if rad 
        Rx = rad2deg(Rx);
        Ry = rad2deg(Ry);
        Rz = rad2deg(Rz);
    end
    
    % Handle singularity?  Consider removing...
    if Rz == 180
        Rz = 179.99;
    end
    if Rx == 180
        Rx = 179.99;
    end
    if Ry == 180
        Ry = 179.99;
    end
          
    % Basic rotation matrices:
    Tx = [1       0           0;
          0  cosd(Rx)  -sind(Rx);
          0  sind(Rx)   cosd(Rx)];

    Ty = [ cosd(Ry) 0  sind(Ry);
               0     1      0;
          -sind(Ry) 0  cosd(Ry)];
      
    Tz = [cosd(Rz)  -sind(Rz)  0;
          sind(Rz)   cosd(Rz)  0;
               0           0   1];

    A = Tz*Ty*Tx;
end