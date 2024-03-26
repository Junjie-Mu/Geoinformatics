%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Positioning and Location Based Services
% A.A. 2023/2024
% Exercise 2:  GPS orbits computation
% 
% Mu Junjie Deng Jianwei Su Jiayi
% 
% References:
%     - http://www.navipedia.net/index.php/Clock_Modelling
%     - http://www.ajgeomatics.com/SatelliteOrbitsRev0828-2010.pdf
%     - http://www.navipedia.net/index.php/GPS_and_Galileo_Satellite_Coordinates_Computation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
set(0,'DefaultFigureWindowStyle','docked');

% Load Almanac of satellite SVN 63, PRN 01 (Block IIR) for 2016
dt0= -7.661711424589D-05;
dt1= -3.183231456205D-12;
dt2=  0.000000000000D+00;
a = 5.153650835037D+03^2;
e = 3.841053112410D-03;
M0 = 1.295004883409D+00;
Omega0 = -2.241692424630D-01;
Omegadot = -8.386063598924D-09;
i0 = 9.634782624741D-01;
idot = -7.286017777600D-11;
w0 = 9.419793734505D-01;
wdot = 0.0;

GMe = 3.986005D+14;
OmegaEdot = 7.2921151467D-05;
% Initialize vector of epochs 
t_0 = 0;     
t_end = 23*3600+59*60+59;
t_step = 30;    
time_epochs = (t_0:t_step:t_end);
% 1) Compute clock offsets and plot it

clock_offset = dt0 + dt1*time_epochs + dt2*time_epochs.^2;
figure(1)
plot(time_epochs, clock_offset, "-");
xlabel('seconds in a day');
ylabel('clock-offset');
title('time-clockoffsets');
ax = gca;
ax.XAxis.Exponent = 0;
% 2) Compute positions in ITRF, X, Y, Z
% Initialize vector of positions
n_epochs = length(time_epochs);
satilte_ORS = zeros(3, n_epochs);
satelite_ITRF = zeros(3, n_epochs);
t=zeros(1, n_epochs);
% compute the mean velocity
n = sqrt(GMe / a^3);

i = 1;
for Dt = t_0 : t_step : (t_end - t_0)
    t(i) = Dt;
    % compute psi
    % a.mean anomaly
    M_t = M0 + n*Dt;
    % b.eccentric anomaly
    E = ecc_anomaly(M_t, e);
    % c.true anomaly
    psi = atan2(sqrt(1-e^2)*sin(E),(cos(E)-e));
    % compute r 
    rt = a*(1-e^2)/(1+e*cos(psi));
    % satellite	coordinates	in	the	ORS
    xt = rt*cos(psi);
    yt = rt*sin(psi);
    satilte_ORS(:,i) = [xt;yt;0];
    % compute w(t)=w0+wdot(t-t0)
    wt = w0 + wdot*Dt;
    % compute i(t)
    it = i0 + idot*Dt;
    % compute W(t)
    omega = Omega0 + (Omegadot-OmegaEdot)*Dt;
    %fill three rotation matrices
    R3o = [cos(-omega),sin(-omega),0;-sin(-omega),cos(-omega),0;0,0,1];
    R1i = [1,0,0;0,cos(-it),-sin(-it);0,sin(-it),cos(-it)];
    R3w = [cos(-wt),sin(-wt),0;-sin(-wt),cos(-wt),0;0,0,1];
    %rotate from OCRS to ITRF
    satelite_ITRF(:,i) = R3o*R1i*R3w*satilte_ORS(:,i);
    i = i+1;
end
% 3) Convert X(t), Y(t), Z(t) to phi(t), lat(t), h(t)
% compute the average radius
Geo = cart2geo(satelite_ITRF);

% 4) Plot satellite's daily trajectory with basemap
figure(2);
% H = subplot(m,n,p), or subplot(mnp), breaks the Figure window
% into an m-by-n matrix of small axes

% Plot groundtracks
% axesm Define map axes and set map properties
ax = axesm ('eqdcylin', 'Frame', 'on', 'Grid', 'on', 'LabelUnits', 'degrees', 'MeridianLabel', 'on', 'ParallelLabel', 'on', 'MLabelParallel', 'south');
% geoshow Display map latitude and longitude data 
%  DISPLAYTYPE can be 'point', 'line', or 'polygon' and defaults to 'line'
geoshow('landareas.shp', 'FaceColor', 'black');
hold on
lon = Geo(1:end,2);
lat = -Geo(1:end,1);
geoshow(lat,lon, 'DisplayType', 'point', 'MarkerEdgeColor', 'green');
% axis EQUAL  sets the aspect ratio so that equal tick mark
% increments on the x-,y- and z-axis are equal in size.
% axis TIGHT sets the axis limits to the range of the data.
axis equal; axis tight;

% Plot height of the satellite 
figure(3)
plot(t(1:10:end), Geo(1:10:end,3), '.g');
%title(['ellipsoidic height variations [km] around mean height = ' num2str(mean(##YourVariables) ' km']);
xlabel('seconds in one day (00:00 - 23:59 = 86400 sec)');
ylabel('[km]');
xlim([1 t(end)]);


% SYNTAX
%  *  [E] = ecc_anomaly(M, e)
% OUTPUT
%  *  E = eccentricity anomaly
% INPUT
%  *  M = mean anomaly
%  *  e = eccentricity of the orbit

function [E] = ecc_anomaly(M, e)
    E = M;
    
    max_iter = 12; %it was 10 when using only GPS (convergence was achieved at 4-6 iterations);
                   % now it set to 12 because QZSS PRN 193 can take 11 iterations to converge
    i = 1;
    dE = 1;
    
    while ((dE > 1e-12) && (i < max_iter))
       E_tmp = E;
       E = M + e * sin(E);
       dE = mod(E - E_tmp, 2 * pi);
       i = i + 1;
    end
    
    if (i == max_iter)
        fprintf('WARNING: Eccentric anomaly does not converge.\n')
    end
    
    E = mod(E, 2 * pi);
end

function result = cart2geo(input)
    n = size(input);
    for i = 1:n(1,2)
        X = input(1,i);
        Y = input(2,i);
        Z = input(3,i);
    
        a = 6378137;
        b = 6356752.314;
        eb = sqrt((a^2-b^2)/b^2);
        r = sqrt(X^2+Y^2);
        e = sqrt(1 - b^2/a^2);
    
        Psi = atan((Z*a)/(r*b));
        Lambda = atan2(Y,X);
        Fai = atan2((Z+eb^2*b*sin(Psi)^3),(r-e^2*a*cos(Psi)^3));
        Rn = a / sqrt(1-e^2*sin(Fai)^2);
        h = (r/cos(Fai) -Rn)/1000 -20189.3326;
    
        result(i, :) = [Fai/pi*180; Lambda/pi*180; h];
    end

end






