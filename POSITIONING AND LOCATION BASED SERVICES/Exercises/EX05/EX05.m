% Positioning & Location Based Services
% A.A. 2023/2034
% EX05: DD analysis and cycle slip identification and repairing
% Author: Mu Junjie
% -------------------------------------------------------------------------
% Guidelines
% 
% Input data = 'CycleSlipsData.txt' text file with
%              col1 = epoch [s]
%              col2 = observed DD [m]
%              col3 = approx DD [m]
%
% 1) Import data in Matlab and graph observed DDs
% 2) For each epoch, compute differences between observed DDs and
%    approximated DDs (residual DDs) and graph them
% 3) For all the couples of consecutive epochs, compute differences between residual DDs 
%    and graph them
% 4) Identify cycle slips and repair them
% 5) Graph the repaired DDs
%
%
% Hints 
%  In step 1) use function 'importdata': 
%  In step 3) try to use function 'diff' instead of a for cycle
%
%  In step 4) use the algorithm explained in the lectures. In
% computation, use 19 cm for wavelenght and 0.20 cycle (3.8 cm) as
% threshold for cycle slips identification and repairing.

clear
close all

% Import data in Matlab and graph observed DDs
[newdata] = importdata('CycleSlipsDataSun.txt', '	', 1);
idata = newdata.data;
disp 'hallo'

threshold = 3.8*1e-2; % [m]
lam = 19*1e-2; % [m]

epoch = idata(:, 1);
observed_DD = idata(:, 2);
approx_DD = idata(:, 3);

figure(1);
plot(epoch, observed_DD, '-');
xlabel('Epoch [s]');
ylabel('Observed DD [m]');
title('Observed DDs');

% For each epoch compute differences between observed DDs and approximated DDs (residual DDs) and graph them
residual_DD = observed_DD - approx_DD;

figure(2);
plot(epoch, residual_DD, '-');
xlabel('Epoch [s]');
ylabel('Residual DD [m]');
title('Residual DDs');

% Compute differences between consecutive epochs of residual DDs (hint: diff or for cycle) and graph them
diff_DD = diff(residual_DD);

figure(3);
plot(epoch(2:end), diff_DD, '-');
xlabel('Epoch [s]');
ylabel('Residual DD Difference [m]');
title('Residual DD Differences');
% Identify cycle slips and repair them (hint: just one for cycle with both the actions)
observed_DD_corr = observed_DD;
for i = 1:length(diff_DD)
    if(abs(diff_DD(i)) > threshold)
        x = diff_DD(i)/lam;
        n = round(x);
        if(lam*abs(n-x)<=threshold)
            observed_DD_corr(i+1:end) = observed_DD(i+1:end) - lam*n;
        end
    end
   
end
% Graph the corrected DDs, the corrected residuals DDs and their differences in time
figure(4)
plot(epoch, observed_DD_corr, '-b','LineWidth', 1)
figure(5);
residual_DD_corr = observed_DD_corr - approx_DD;
plot(epoch, residual_DD_corr, '-b','LineWidth', 1)
figure(6);
diff_DD_corr = diff(residual_DD_corr);
plot(epoch(2:end), diff_DD_corr, '-g','LineWidth', 0.7)
xlabel('Epoch [s]');
legend('corrected residuals DDs','differences')