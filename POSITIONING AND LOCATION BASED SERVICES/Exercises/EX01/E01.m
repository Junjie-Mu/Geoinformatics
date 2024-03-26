%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Positioning and Location Based Services
% A.A. 2023/2024
% Exercise 1:  Reference Frames
% 
% Mu Junjie Deng Jianwei Su Jiayi
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all

%1 Import Data
Lat_O = [44, 23, 24];
Lon_O = [8, 56, 20];
H_O = 70;

%2 Origin ITRF to Global Cartesian
a = 6378137;
f = 1/298.257222100882711243;
e = sqrt(f*(2 - f));
Rn = a / sqrt(1-e^2*sin(ToRad(Lat_O))^2);
Lambda0 = ToRad(Lon_O);
Fai0 = ToRad(Lat_O);
X0_GC = (Rn + H_O)*cos(Fai0)*cos(Lambda0);
Y0_GC = (Rn + H_O)*cos(Fai0)*sin(Lambda0);
Z0_GC = (Rn*(1-e^2)+H_O)*sin(Fai0);
P0 = [X0_GC;Y0_GC;Z0_GC];

%3 LL TO LC 
qsi = ToRad([0,0,10.23]);
eta =ToRad([0,0,9.5]);
alpha = ToRad([30, 27, 18]);
Rx_qsi = [1,0,0; 0,cos(qsi),-sin(qsi); 0,sin(qsi),cos(qsi)];
Ry_eta = [cos(eta),0,-sin(eta); 0,1,0; sin(eta),0,cos(eta)];
Rz_alpha = [cos(alpha),sin(alpha),0; -sin(alpha),cos(alpha),0; 0,0,1];

RLC2LL = Rz_alpha*Ry_eta*Rx_qsi;
RLL2LC =transpose(RLC2LL);
A = [0;30;0];
B = [0;-30;0];
C = [200;0;0];

A_LC = RLL2LC * A;
B_LC = RLL2LC * B;
C_LC = RLL2LC * C;

%4 LC To ITRF GC
R0 = [-sin(Lambda0), cos(Lambda0), 0;
      -sin(Fai0)*cos(Lambda0), -sin(Fai0)*sin(Lambda0), cos(Fai0);
      cos(Fai0)*cos(Lambda0), cos(Fai0)*sin(Lambda0), sin(Fai0)];
A_GC = P0 + transpose(R0)*A_LC;
B_GC = P0 + transpose(R0)*B_LC;
C_GC = P0 + transpose(R0)*C_LC;

%5 ITRF GC to ETRF GC throuth EPN
ETRF_A = [4509854.8133, 709344.7333, 4439228.7611];
ETRF_B = [4509885.8305, 709380.3976, 4439191.8018];
ETRF_C = [4509773.4717, 709521.8569, 4439282.7125];

%6 ETRF GC to Geodatic
A_Geo = ToGeo(ETRF_A);
B_Geo = ToGeo(ETRF_B);
C_Geo = ToGeo(ETRF_C);

%7 Accuracies from LL to LC
deltaA = 0.02;
deltaC = 0.1;
CGCO=[0.1^2,0,0;0,0.1^2,0;0,0,0.1^2];
SD_AB = GetConvariance(deltaA, RLL2LC,R0);
SD_C = GetConvariance(deltaC, RLL2LC,R0);

%8 Save on txt
fileID = fopen('result','w');

fprintf(fileID, 'Group Member: Mu Junjie, Deng Jianwei, Su Jiayi\n');
fprintf(fileID, 'Local cartesian coordinates of points A, B, C in meters\n');
fprintf(fileID, 'Point A = \nE:  %.3f\nN:  %.3f\nU:  %.3f\n',A_LC(1),A_LC(2),A_LC(3));
fprintf(fileID, 'Point B = \nE:  %.3f\nN:  %.3f\nU:  %.3f\n',B_LC(1),B_LC(2),B_LC(3));
fprintf(fileID, 'Point C = \nE:  %.3f\nN:  %.3f\nU:  %.3f\n',C_LC(1),C_LC(2),C_LC(3));

fprintf(fileID, '\nITRF global cartesian coordinates of points A,B,C \n');
fprintf(fileID, 'Point A = \nX:  %.3f\nY:  %.3f\nZ:  %.9f\n',A_GC(1),A_GC(2),A_GC(3));
fprintf(fileID, 'Point B = \nX:  %.3f\nY:  %.3f\nZ:  %.9f\n',B_GC(1),B_GC(2),B_GC(3));
fprintf(fileID, 'Point C = \nX:  %.3f\nY:  %.3f\nZ:  %.9f\n',C_GC(1),C_GC(2),C_GC(3));

fprintf(fileID, '\nETRF cartesian coordinates of points A, B, C  \n');
fprintf(fileID, 'Point A = \nX:  %.4f\nY:  %.4f\nZ:  %.4f\n',ETRF_A(1),ETRF_A(2),ETRF_A(3));
fprintf(fileID, 'Point B = \nX:  %.4f\nY:  %.4f\nZ:  %.4f\n',ETRF_B(1),ETRF_B(2),ETRF_B(3));
fprintf(fileID, 'Point C = \nX:  %.4f\nY:  %.4f\nZ:  %.4f\n',ETRF_C(1),ETRF_C(2),ETRF_C(3));

fprintf(fileID, '\nETRF geodetic coordinates of points A,B,C   \n');
fprintf(fileID, 'Point A = \nlat:  [%d, %d, %.4f]\nlon:  [%d, %d, %.4f]\nh:    [%.3f]\n',A_Geo(1,:),A_Geo(2,:),A_Geo(3,3));
fprintf(fileID, 'Point B = \nlat:  [%d, %d, %.4f]\nlon:  [%d, %d, %.4f]\nh:    [%.3f]\n',B_Geo(1,:),B_Geo(2,:),B_Geo(3,3));
fprintf(fileID, 'Point C = \nlat:  [%d, %d, %.4f]\nlon:  [%d, %d, %.4f]\nh:    [%.3f]\n',C_Geo(1,:),C_Geo(2,:),C_Geo(3,3));

fprintf(fileID, '\nStandard deviations of points A,B in East, North, Up in cm   \n');
fprintf(fileID, 'E:  %.1f\nN:  %.1f\nU:  %.1f\n',SD_AB(1),SD_AB(2),SD_AB(3));
fprintf(fileID, '\nStandard deviations of points c in East, North, Up in cm   \n');
fprintf(fileID, 'E:  %.1f\nN:  %.1f\nU:  %.1f\n',SD_C(1),SD_C(2),SD_C(3));

%Functions
%Convert to radian
function result = ToRad(input)
    a =input(1);
    b =input(2);
    c =input(3);
    result = (a+b/60+c/3600)*pi/180;
end

%Cartesian to Geodetic
function result = ToGeo(input)
    X = input(1);
    Y = input(2);
    Z = input(3);

    a = 6378137;
    b = 6356752.314;
    eb = sqrt((a^2-b^2)/b^2);
    r = sqrt(X^2+Y^2);
    e = sqrt(1 - b^2/a^2);

    Psi = atan((Z*a)/(r*b));
    Lambda = atan2(Y,X);
    Fai = atan2((Z+eb^2*b*sin(Psi)^3),(r-e^2*a*cos(Psi)^3));
    Rn = a / sqrt(1-e^2*sin(Fai)^2);
    h = r/cos(Fai) -Rn;

    result = [ToSexagesimal(Fai/pi*180);
              ToSexagesimal(Lambda/pi*180);
              [0,0,h]];
end

%Convert to degree
function result = ToSexagesimal(input)
    degreesInt = fix(input);  
    minutesFloat = (input - degreesInt) * 60;  
    minutesInt = fix(minutesFloat); 
    secondsFloat = (minutesFloat - minutesInt) * 60;  
    result = [degreesInt, minutesInt, secondsFloat];
end

%Convariance of GC
function result = GetConvariance(delta,RLL2LC,R0)
    C_GC_O = diag([0.1^2, 0.1^2, 0.1^2]);
    C_LL = diag([delta^2,delta^2,delta^2]);
    C_LC = RLL2LC*C_LL*transpose(RLL2LC);
    C_GC = C_GC_O + transpose(R0)*C_LC*R0;
    C = R0*C_GC*transpose(R0);
    result = 100*[sqrt(C(1,1)),sqrt(C(2,2)),sqrt(C(3,3))];
end