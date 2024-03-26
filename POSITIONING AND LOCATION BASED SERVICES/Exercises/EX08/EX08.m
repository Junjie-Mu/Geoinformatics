%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EX08: Inertial Navigation
% POS&LBS A.A. 2023/2024
% MU JUNJIE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc;
% read database of control points and relevant vector of measurements
control_points = importdata("control_points_db.txt");
% read user points and relevant measurements
user_points = importdata("user_db.txt");
% idem

% for each user point
for i = 1:height(user_points)
    % for each control point
    for j = 1:height(control_points)
            % find the distance of measurements of the user from the
            % measurements of the control point: (m_u-m_cp)*(m_u-m_cp)
            diff(j,:) = user_points(i,2:end) - control_points(j,4:end);
    end
    vect = sum(diff.^2,2);
    % identify the minimum
    cp_id = find(vect == min(vect));
    % attribute to the user the position of the control point with the
    % nearest measurements
    nearest_points(i,:) = [i, control_points(cp_id,1:3)];
end
% some nice plot with the movements of the user
figure;
rectangle('Position', [-5, -5, 30, 20]);
hold on;
rectangle('Position', [0, 0, 20, 10]);
hold on;
plot(nearest_points(:,3), nearest_points(:,4), '-g','DisplayName','user');
hold on;
for x = 2:2:18
    for y = 2:2:8
        plot(x, y, 'ko', 'MarkerSize', 1.5);
    end
end
plot(nearest_points(:,3), nearest_points(:,4), '-g','DisplayName','user');
legend('user','control points')
hold off;
title('User Trajectory');
