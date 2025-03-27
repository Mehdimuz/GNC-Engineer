% By Benjamin Reifler
% Implements the IAU-2000/2006 reduction, partially adapted from MATLAB aero toolbox
% Inputs:   JD (Julian date)
%           coefs (1600x17,1275x17,66x11)
% Outputs:  parameters

function [W,Q] = eci2ecef_INIT(JD,coefs)
    % convert JD to modified JD
    JD = JD - 2400000.5;
    jd = JD - 51544.5;

    % precompute powers of t (Julian centuries since epoch)
    t = jd/36525;
    t2 = t*t;
    t3 = t2*t;
    t4 = t3*t;
    t5 = t4*t;

    % etc.
    mMoon = 485868.249036 + 1717915923.2178*t + 31.8792*t2 + 0.051635*t3 - 0.00024470*t4;
    mSun = 1287104.793048 + 129596581.0481*t - 0.5532*t2 + 0.000136*t3 - 0.00001149*t4;
    umMoon = 335779.526232 + 1739527262.8478*t - 12.7512*t2 - 0.001037*t3 + 0.00000417*t4;
    dSun = 1072260.703692 + 1602961601.2090*t - 6.3706*t2 + 0.006593*t3 - 0.00003169*t4;
    omegaMoon = 450160.398036 - 6962890.5431*t + 7.4722*t2 + 0.007702*t3 - 0.00005939*t4;

    lMercury = 4.402608842 + 2608.7903141574*t;
    lVenus = 3.176146697 + 1021.3285546211*t;
    lEarth = 1.753470314 + 628.3075849991*t;
    lMars = 6.203480913 + 334.06124267*t;
    lJupiter = 0.599546497 + 52.9690962641*t;
    lSaturn = 0.874016757 + 21.329910496*t;
    lUranus = 5.481293872 + 7.4781598567*t;
    lNeptune = 5.311886287 + 3.8133035638*t;
    pa = 0.02438175*t + 0.00000538691*t2;

    nutationV = mod([[mMoon; mSun; umMoon; dSun; omegaMoon]/3600*pi/180; lMercury; lVenus; lEarth; lMars; lJupiter; lSaturn; lUranus; lNeptune; pa],2*pi);

    X0 = -16617 + 2004191898*t - 429782.9*t2 - 198618.34*t3 + 7.578*t4 + 5.9285*t5;
    Y0 = -6951 - 25896*t - 22407274.7*t2 + 1900.59*t3 + 1112.526*t4 + 0.1358*t5;

    S0 = 94 + 3808.65*t - 122.68*t2 - 72574.11*t3 + 27.98*t4 + 15.62*t5;

    FX = zeros(length(coefs.X),1);
    FX(1:1306,:) = ones(1306,1);
    FX(1307:1559,:) = repmat(t,[253 1]);
    FX(1560:1595,:) = repmat(t2,[36 1]);
    FX(1596:1599,:) = repmat(t3,[4 1]);
    FX(1600,:) = t4;
    argX = coefs.X(:,4:17)*nutationV;
    X = sum((coefs.X(:,2).*sin(argX(:,1)) + coefs.X(:,3).*cos(argX(:,1))).*FX(:,1));

    FY = zeros(length(coefs.Y),1);
    FY(1:962,:) = ones(962,1);
    FY(963:1239,:) = repmat(t,[277 1]);
    FY(1240:1269,:) = repmat(t2,[30 1]);
    FY(1270:1274,:) = repmat(t3,[5 1]);
    FY(1275,:) = t4;
    argY = coefs.Y(:,4:17)*nutationV;
    Y = sum((coefs.Y(:,2).*sin(argY(:,1)) + coefs.Y(:,3).*cos(argY(:,1))).*FY(:,1));

    FS = zeros(length(coefs.S),1);
    FS(1:33,:) = ones(33,1);
    FS(34:36,:) = repmat(t,[3 1]);
    FS(37:61,:) = repmat(t2,[25 1]);
    FS(62:65,:) = repmat(t3,[4 1]);
    FS(66,:) = t4;
    argS = coefs.S(:,4:11)*[nutationV(1:5,:);nutationV(7:8,:);nutationV(14,:)];
    S = sum((coefs.S(:,2).*sin(argS(:,1)) + coefs.S(:,3).*cos(argS(:,1))).*FS(:,1));

    X = X + X0;
    Y = Y + Y0;
    S = S + S0;

    X = X*1e-6/3600*pi/180;
    Y = Y*1e-6/3600*pi/180;

    S = S*1e-6/3600*pi/180 - X*Y/2;
    E = atan2(Y,X);
    d = atan(sqrt((X.^2+Y.^2)./(1-X.^2-Y.^2)));

    % Exported Nutation Values:
    cQ  = cos(E);
    sQ  = sin(E);
    cQ2 = cos(d);
    sQ2 = sin(d);
    cQ3 = cos(-E-S);
    sQ3 = sin(-E-S);

    % Exported Polar Motion Values:
    sp = -0.000047*t/3600*pi/180;
    cW = cos(sp);
    sW = sin(sp);

    % polar motion and nutation matrices:
    W = [cW sW 0 ; -sW cW 0 ; 0 0 1];
    Q = [cQ3 sQ3 0 ; -sQ3 cQ3 0 ; 0 0 1]*[cQ2 0 -sQ2 ; 0 1 0 ; sQ2 0 cQ2]*[cQ sQ 0 ; -sQ cQ 0 ; 0 0 1];
end
