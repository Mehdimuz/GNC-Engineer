function [euler] = a2e(A,rad)  
    if nargin == 1
        rad = false;
    end
    Ry = -asin(A(3,1));
    Rx = atan2(A(3,2), A(3,3));
    Rz = atan2(A(2,1), A(1,1));
    
    % TODO: Fix this (low priority)
%     if Ry == 0
%         Rx  = atan2(-A(2,1), A(1,1));
%     elseif Ry == pi
%         Rx = atan2(A(2,1), A(1,1));
%     end
    
%     if Rx < 0
%         Rx = Rx + 2*pi; 
%     end
%     if Rz < 0
%         Rz = Rz + 2*pi;
%     end
    
    if rad
        euler = [Rx; Ry; Rz];
    else
        euler = [rad2deg(Rx); rad2deg(Ry); rad2deg(Rz)];
    end
end