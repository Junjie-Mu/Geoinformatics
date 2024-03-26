% Positioning & Location Based Services
% A.A. 2023/2034
% EX06
% Author: Mu Junjie
%

clear
close all

%Import data from RTK
%Next data replace by s4_t1_2.csv
RTK_data1 = readtable("s4_t1.csv");
data1_northing = RTK_data1.Northing;
data1_easting = RTK_data1.Easting;
n = height(RTK_data1);
utmzone = repmat('32 N', n, 1);
%Transfer UTM to degree
[Lat, Lon] = utm2deg(data1_easting,data1_northing,utmzone);
%Transfer UTK dgree to csv
data1_table = table(Lat, Lon);
csv_path = 'UTK_data.csv';
writetable(data1_table, csv_path);

%Import data from Matlab
%Next data replace by S4-T1-2.mat
mat_data = load("S4-T1.mat");
mat_latitude = mat_data.Position.latitude;
mat_longitude = mat_data.Position.longitude;
%Transfer mat data to csv
data_table = table(mat_latitude, mat_longitude);
csv_path = 'mat_data.csv';
writetable(data_table, csv_path);
[E_mat,N_mat,utm_zone] = deg2utm(mat_latitude,mat_longitude);

%Plot the data on the webmap using geocoordinates
webmap("worldstreetmap")
geoplot(Lat,Lon,'bo-','MarkerSize',5,'LineWidth',2);
hold on 
geoplot(mat_latitude,mat_longitude,'bo-','Color','red','MarkerSize',5,'LineWidth',2);
hold off
%Add a title
title('Stonex vs Matlib data');

%Estimate errors
distance = zeros(size(E_mat));
for i = 1:length(E_mat)
    dist_allPoints = sqrt((data1_easting-E_mat(i)).^2 + (data1_northing-N_mat(i)).^2);
    distance(i) = min(dist_allPoints);
end

%Plots of errors
figure(2)
plot((1:length(mat_latitude)),distance)

%Print the statistics(mean,max,min STD)
mean = mean(distance);
max = max(distance);
min = min(distance);
STD = std(distance);

tab = table(mean,max,min,STD);
writetable(tab,'statistics.csv')
