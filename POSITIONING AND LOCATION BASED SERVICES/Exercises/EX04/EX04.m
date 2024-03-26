% Positioning and Location Based Services
% A.A. 2023/2024
% 4rd Exercise
%
% @author: Mu Junjie
clc
clear all

%getting obs values
run('Ex04_variables.m')
                                                                                               
%initialization of parameters
th_convergence = 0.1;
%max number of iterations
max_iter = 10;
%threshold for convergence
X_init = [0 0 0 0]; %approx coordinates of the receiver and clock offset
% storage of iterate results
Xr_est = zeros(4);
% Least squares iteration
Xr = X_init;
for i=1:max_iter
    [lat, lon, h] = cart2geo(Xr(1,1),Xr(1,2),Xr(1,3));
    % topocentric positions of all satellites
    n_sat = size(xyz_sat,1);
    %LS design matrix
    A = zeros(n_sat,4);
    %LS known term vector
    b = zeros(n_sat,1);
    for j = 1:n_sat
        Xs = xyz_sat(j,:);
        [Az, El, D] = topocent(Xr,Xs);
        delta_xyz  = Xs - Xr(1:3);
        % approximate geodetic coordinates of the receiver
        
        % tropospheric and ionospheric corrections
        iono_correct = iono_correction(lat,lon,Az,El,time_rx,ionoparams);
        tropo_correc = tropo_correction(h,El);
        % LS known term
        b(j) = iono_correct + tropo_correc + D - s_light*dtS(j);
        % LS A matrix
        A(j,1:3) = -delta_xyz/D;
        A(j,4) = 1;
    end
    % L+east square solution for the corrections to the apriori
    y0 = pr_C1-b;
    N = transpose(A)*A;
    inv_N = inv(N);
    % Estimated coordinates of the receiver: 
    X0 = inv_N*transpose(A)*y0;
    Xr_est = transpose(X0);
    % approximate + estimated correction
    Xr_est(1,1:3) = Xr(1,1:3) + Xr_est(1:1:3);
    Xr_pre = Xr;
    Xr = Xr_est;
    delta = Xr - Xr_pre;
    %check convergence of the result and, in case exit
    if max(abs(delta(1:3))) < th_convergence
        break
    end
    % check at the end that convergence did not fail
    if i == max_iter
        disp('Convergence failed');
    end
end

% final estimated unknowns
R_pos = Xr_est;
R_pos(1, 4) = R_pos(1, 4) / s_light;
% LS residuals and sigma2
residuals = pr_C1 - A * transpose(R_pos) - b;
sigma2 = transpose(residuals)*residuals / (n_sat - 4);
% covariance matrix of the estimated coordinates
Cxx = sigma2 * inv(N);
Qxx = transpose(inv_N(1:3, 1:3));
% Rotate and PDOP
[Fai0, Lambda0, h0] = cart2geo(Xr(1,1),Xr(1,2),Xr(1,3));
R0 = [-sin(Lambda0), cos(Lambda0), 0;
      -sin(Fai0)*cos(Lambda0), -sin(Fai0)*sin(Lambda0), cos(Fai0);
      cos(Fai0)*cos(Lambda0), cos(Fai0)*sin(Lambda0), sin(Fai0)];
Qll = R0*Qxx*transpose(R0);
PDOP = sqrt(Qll(1,1)+Qll(2,2)+Qll(3,3));
% print results
i_print = sprintf('The total number of iterations is: %d', i);
disp(i_print);
coord_print = sprintf('The coordinates of the receiver are: [%d %d %d]', R_pos(1,1:3));
disp(coord_print);
offset_print = sprintf('The clock offset of the receiver is: %d', R_pos(1, 4));
disp(offset_print);
PDOP_print = sprintf('PDOP value is %f', PDOP);
disp(PDOP_print);

%%

%% repeat with CutOfAngle

CutOfAngle = 5;

xyz_sat_co = zeros(0, 3);
pr_C1_co = zeros(0); 
dtS_co = zeros(0);

% extract satellites above cut off
for j=1:n_sat
    Xs = xyz_sat(j,:);
    [Az, El_co, D] = topocent(Xr,Xs);
    if abs(El_co) > CutOfAngle
        % store coordinates, otherwise exclude
        xyz_sat_co = [xyz_sat_co;xyz_sat(j, :)];
        pr_C1_co = [pr_C1_co;pr_C1(j)];
        dtS_co = [dtS_co;dtS(j)];
    end
end

% repeat same computations
for i=1:max_iter
    [lat, lon, h] = cart2geo(Xr(1,1),Xr(1,2),Xr(1,3));
    % topocentric positions of all satellites
    n_sat = size(xyz_sat_co,1);
    %LS design matrix
    A = zeros(n_sat,4);
    %LS known term vector
    b = zeros(n_sat,1);
    for j = 1:n_sat
        Xs = xyz_sat_co(j,:);
        [Az, El, D] = topocent(Xr,Xs);
        delta_xyz  = Xs - Xr(1:3);
        % approximate geodetic coordinates of the receiver
        
        % tropospheric and ionospheric corrections
        iono_correct = iono_correction(lat,lon,Az,El,time_rx,ionoparams);
        tropo_correc = tropo_correction(h,El);
        % LS known term
        b(j) = iono_correct + tropo_correc + D - s_light*dtS_co(j);
        % LS A matrix
        A(j,1:3) = -delta_xyz/D;
        A(j,4) = 1;
    end
    % L+east square solution for the corrections to the apriori
    y0 = pr_C1_co-b;
    N = transpose(A)*A;
    inv_N = inv(N);
    % Estimated coordinates of the receiver: 
    X0 = inv_N*transpose(A)*y0;
    Xr_est = transpose(X0);
    % approximate + estimated correction
    Xr_est(1,1:3) = Xr(1,1:3) + Xr_est(1:1:3);
    Xr_pre = Xr;
    Xr = Xr_est;
    delta = Xr - Xr_pre;
    %check convergence of the result and, in case exit
    if max(abs(delta(1:3))) < th_convergence
        break
    end
    % check at the end that convergence did not fail
    if i == max_iter
        disp('Convergence failed');
    end
end

% final estimated unknowns
R_pos_co = Xr_est;
R_pos_co(1, 4) = R_pos_co(1, 4) / s_light;
% LS residuals and sigma2
residuals = pr_C1_co - A * transpose(R_pos_co) - b;
sigma2 = transpose(residuals)*residuals / (n_sat - 4);
% covariance matrix of the estimated coordinates
Cxx = sigma2 * inv(N);
Qxx = transpose(inv_N(1:3, 1:3));
% Rotate and PDOP
[Fai0, Lambda0, h0] = cart2geo(Xr(1,1),Xr(1,2),Xr(1,3));
R0 = [-sin(Lambda0), cos(Lambda0), 0;
      -sin(Fai0)*cos(Lambda0), -sin(Fai0)*sin(Lambda0), cos(Fai0);
      cos(Fai0)*cos(Lambda0), cos(Fai0)*sin(Lambda0), sin(Fai0)];
Qll = R0*Qxx*transpose(R0);
PDOP_co = sqrt(Qll(1,1)+Qll(2,2)+Qll(3,3));

% print
disp('-------------------------------------------------------------');
COA_print = sprintf('Cut off Angle: %d',CutOfAngle);
disp(COA_print);
i_print = sprintf('The total number of iterations (COA) is: %d', i);
disp(i_print);
coord_print = sprintf('The coordinates of the receiver (COA) are: [%d %d %d]', R_pos_co(1,1:3));
disp(coord_print);
offset_print = sprintf('The clock offset of the receiver (COA) is: %d', R_pos_co(1,4));
disp(offset_print);
PDOP_print = sprintf('PDOP (COA) value is %f', PDOP_co);
disp(PDOP_print);