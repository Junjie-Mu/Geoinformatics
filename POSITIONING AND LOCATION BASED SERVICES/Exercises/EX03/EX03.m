% Positioning and Location Based Services
% A.A. 2023/2024
% 3rd Exercise: Ionospheric delay computation
%
% @author: Mu Junjie

%load parameters
line1 = 'GPSA 7.4506D-09 1.4901D-08 -5.9605D-08 -1.1921D-07 IONOSPHERIC CORR';
line2 = 'GPSB 9.2160D+04 1.3107D+05 -6.5536D+04 -5.2429D+05 IONOSPHERIC CORR';
ionoparams = [cell2mat(textscan(line1, '%*s %f %f %f %f %*s')) ...
cell2mat(textscan(line2, '%*s %f %f %f %f %*s'))];
%initialize values for the zenital cycle
az1 = 0;
el1 = 90;
%initialize matrix
Iono_map1 = zeros();
%time, phi and lambda cycle
time1 = [0, 6*3600, 12*3600, 18*3600];
phi1 = (-80:0.5:80);
lambda1 = (-180:0.5:180);
% plots
%[phi_grid,lambda_grid] = meshgrid(phi1,lambda1);
[lambda_grid,phi_grid] = meshgrid(lambda1,phi1);
figure(1)
title('Ionospheric Error Maps')
for i = 1:length(time1)
    for s  = 1:length(phi1)
        for t = 1:length(lambda1)
         Iono_map1(s,t,i) = iono_correction(phi1(s),lambda1(t),az1,el1,time1(i),ionoparams);
        end
    end
    subplot(2,2,i);
    Iono_map0 = Iono_map1(:,:,1);
    geoshow(phi_grid, lambda_grid, Iono_map1(:,:,i), 'DisplayType','texturemap','facealpha',.5)
    hold on
    geoshow('landareas.shp', 'FaceColor', 'none');
    title(['time = ', num2str(time1(i)/3600),':00']);
    xlabel('longitude [deg]')
    ylabel('latitude [deg]')
    xlim([-180 180]);
    ylim([-80 80]);
    colormap(jet);
end
hp4 = get(subplot(2,2,4),'Position');
colorbar('Position', [hp4(1)+hp4(3)+0.028  hp4(2)  0.03  hp4(2)+hp4(3)*2.1]);

%%polar map in Milano
% Milano position in degrees
phi2 = 45 + 28 / 60 + 38.28 / 60^2; %degrees
lambda2 = 9 + 10 / 60 + 53.40 / 60^2; %degrees
%inizialize values for the cycle

% matrix inizialization
Iono_map2 = zeros();
%time, elevation and azimuth cycle 
time2 = [0,12*3600];
azimuth2 = (-180:0.5:180);
elevation2 = (0:0.5:90);
%plots
[Az, El] = meshgrid(azimuth2, elevation2);
for i = 1 : length(time2)
    figure(i + 1)
    title(['Ionospheric Error Polar Map for Milan Observer time = ', num2str(time2(i)/3600),':00']);
    axesm('eqaazim', 'MapLatLimit', [0 90]);
    hold on
    axis off
    framem on
    gridm on

    mlabel on
    plabel on;
    setm(gca,'MLabelParallel',0)
    for s  = 1:length(azimuth2)
        for t = 1:length(elevation2)
         Iono_map2(t,s,i) = iono_correction(phi2,lambda2,azimuth2(s),elevation2(t),time2(i),ionoparams);
        end
    end
    geoshow(El, Az, Iono_map2(:,:,i), 'DisplayType','texturemap', 'facealpha',.6)
    colormap(jet)

    hcb = colorbar('eastoutside');
    set(get(hcb,'Xlabel'),'String','Legend')
end
