% Positioning & Location Based Services
% A.A. 2023/2034
% EX07
% Author: Mu Junjie
%

clc
clear all
%load data without errors
data = load("Inertial_data.dat");
epoch = data(:, 1);
%accelerations
acc = data(:, 2:3); 
%omega in zed
omegaz = data(:,4);

%load data with errors (repeat)
data_ni = load("Inertial_data_ni.dat");
epoch_ni = data_ni(:, 1);
%accelerations
acc_ni = data_ni(:, 2:3); 
%omega in zed
omegaz_ni = data_ni(:,4);

%call function to compute trajectory without errors
traject = CalcTrajectory(acc, omegaz, epoch);
%call function to compute trajectory with errors
traject1 = CalcTrajectory(acc_ni, omegaz_ni, epoch_ni);
%plot comparison on the same plot
figure(1)
plot(traject(1,:), traject(2,:), 'r-', 'LineWidth', 3);
hold on

plot(traject1(1,:), traject1(2,:), 'b-', 'LineWidth', 1);
xlabel('X-IRS');
ylabel('Y-IRS');

function traject = CalcTrajectory(acc, omegaz, epoch)
    %Initialize
    n = size(epoch);
    ax = acc(:,1);
    ay = acc(:,2);
    vx = zeros(n);
    vy = zeros(n);
    dx = zeros(n);
    dy = zeros(n);
    dt = diff(epoch);
    alpha = zeros(n);
    ay_clean = zeros(n);
    XY_IRS = zeros(2,125);
    XY_IRS(:,1) = [100;100];
    for t = 2:length(epoch)
        %Compute X velocities and delta positions in body frame
        vx(t) = vx(t-1) + ax(t)*dt(t-1);
        dx(t) = vx(t)*dt(t-1) + 1/2*ax(t)*dt(t-1)^2;
        %Clean apparent centrifugal from Y accelerate
        ay_clean(t) = ay(t) - omegaz(t)*vx(t);
        %Compute skidding velocity and displacement
        vy(t) = vy(t-1) + ay_clean(t)*dt(t-1);
        dy(t) = vy(t)*dt(t-1) + 1/2*ay_clean(t)*dt(t-1)^2;
        %Compute asset angles alpha
        alpha(t) = alpha(t-1) + omegaz(t)*dt(t-1);
        R_alpha = [cos(alpha(t)),sin(alpha(t));
                   -sin(alpha(t)), cos(alpha(t))];
        dxy = [dx(t);dy(t)];
        XY_IRS(:,t) = XY_IRS(:,t-1) + R_alpha*dxy;
        traject = XY_IRS;
    end
end
