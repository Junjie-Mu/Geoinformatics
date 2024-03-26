% Compute the coordinates of a GPS receiver using one epoch of pseudorange
% observations and the given position of the satellites in view
%
% Meaning of the variables 
%  - time_rx        epoch of observation
%  - pr_C1          array pseudo range observations 
%  - sat_ids        array with the list of the satellites IDs
% 
%  - xyz_sat        estimated positions of the satellites
%  - dtS            estimated clock error of the satellites
%  - ionoparams     8 iono parameters
%
%  - xyz_real       "real" position of the receiver, to check the quality
%                   of your computation
%
% Useful functions:
%  - [az, el, dist] = topocent(xyz_approx, xyz_sat);
%  - [phiR, lamR, hR] = cart2geod(xyz_approx(1), xyz_approx(2), xyz_approx(3));
%  - err_tropo = tropo_error_correction(hR, el);
%  - err_iono = iono_error_correction(phiR, lamR, az, el, time_rx, ionoparams);
%
% Use 10 iterations of linearized LS
% Hint: to compute the satellite ionospheric and troposferic errors
% you need at least a rough estimation of the receiver position.
% Start your iterations using the center of the Earth as approximate
% positions for the receiver.
%
% Compute PDOP, the coordinates of the receiver, the error of the
% estimation, and the covariance matrix of the error
%
% Bonus, try to ignore the satellite below 5 degree of elevation
%

s_light = 299792458; % speed of light

gps_week = 1906;
gps_week_start = 17;
year = 2016;
month = 7;
day = 22;
hour = 2;
minute = 0;
second = 0;

% Convert it in GPS double format
[time_rx] = (day - gps_week_start)*24*3600 + hour*3600 + second*60;

% satellite list
sat_ids = [ 'G28';
            'G05';
            'G13';
            'G07';
            'G20';
            'G09';
            'G08';
            'G02';
            'G21';
            'G30';
            'G15'];

% observed pseudorange of code 
pr_C1 = [22578312.093;   
         20984179.054; 
         22340643.025;
         21815745.186;
         23719962.313;
         24558115.868;
         25751171.706;
         24359848.780;
         26560055.854;
         20547846.341;
         25187331.212];

% satellite positions and clock error at the observation epoch:
% from the function [xyz_sat dtS] = get_sat_pos(obs_head(1:25));
xyz_sat =    1.0e+07 * [
   2.266904303720417   1.376019580610336   0.242665876089085;
   1.793483097238855  -0.684611592059353   1.836848309959032;
   1.237349240389060  -1.488069273671674   1.810628493493055;
   0.682633532041234   1.366381869196001   2.184884029980142;
   0.141020153293916  -1.610493243878792   2.092127478822572;
   0.710758426920100   2.494566375976196   0.565262210580487;
  -0.670387964847087   2.192133222131345   1.342746581584370;
   2.183170669948725  -1.415437238094089  -0.371984191760939;
  -1.019768755765267  -1.243833666228832   2.189467478141541;
   1.528675973015969   0.745824912302640   2.042145267368744;
   0.467596411393501  -2.316970109165663   1.162832980243857];

dtS = [
      0.000538885950029137
     -0.000103714172891042
     -3.26664571204891e-05
      0.000440397108129438
      0.000425625330509237
      0.000171981683578018
     -4.36651382082638e-05
      0.000573964626877986
     -0.000528855944540131
      0.000141099019219313
     -0.000320324134333714];
    
 line1 = 'GPSA   0.7451D-08  0.1490D-07 -0.5960D-07 -0.1192D-06       IONOSPHERIC CORR';
 line2 = 'GPSB   0.9216D+05  0.1311D+06 -0.6554D+05 -0.5243D+06       IONOSPHERIC CORR';
 ionoparams = [cell2mat(textscan(line1, '%*s %f %f %f %f %*s')) ...
     cell2mat(textscan(line2, '%*s %f %f %f %f %*s'))]; 

% Marianna Alghisi
% END OF TEXT -------------------------------------------------------------
