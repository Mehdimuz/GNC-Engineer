% A = dcmeci2ecef_codegen(JD,coefs)
% By Benjamin Reifler
% Implements the IAU-2000/2006 reduction, partially adapted from MATLAB aero toolbox
% Inputs:   JD (Julian date)
%           params
% Outputs:  A (DCM)

function A = eci2ecef(JD,W,Q) %#codegen
    % convert JD to modified JD
    JD = JD - 2400000.5;
    jd = JD - 51544.5;
    jdf = mod(jd, 1);

    % Earth rotation angle
    thetaERA = mod(2*pi*(jdf + 0.7790572732640 + 0.00273781191135448*jd),2*pi);
    cR = cos(thetaERA);
    sR = sin(thetaERA);
    R = [cR sR 0 ; -sR cR 0 ; 0 0 1];
    
    % Full Attitude Matrix:
    A = W(:,:,1)*R(:,:,1)*Q(:,:,1);
end
