function [delay] = iono_correction(lat, lon, az, el, time, ionoparams)

%--------------------------------------------------------------------------
%   KLOBUCHAR MODEL
%   Computation of the pseudorange correction due to ionospheric effect.
%
%   Parameters
%   ----------
%   lat : latitude in degrees.
%   lon : longitude in degrees.
%   az : azimuth in degrees.
%   el : elevation in degrees.
%   GPStime : time in seconds with respect the beginning of the GPS week.
%   ionoparams : vector containing the alpha and bete ionosperic correction 
%                parameters from NAV files.
% 
%   Returns
%   -------
%   T_iono : ionosperic effect in meters.
%--------------------------------------------------------------------------


    %initialization
    delay = zeros(size(el));
    c = 299792458;
    
    %ionospheric parameters
    a0 = ionoparams(1);
    a1 = ionoparams(2);
    a2 = ionoparams(3);
    a3 = ionoparams(4);
    b0 = ionoparams(5);
    b1 = ionoparams(6);
    b2 = ionoparams(7);
    b3 = ionoparams(8);
    
    %elevation from 0 to 90 degrees
    el = abs(el);
    
    %conversion to semicircles
    lat = lat / 180;
    lon = lon / 180;
    az = az / 180;
    el = el / 180;
    
    psi = (0.0137 ./ (el+0.11)) - 0.022;
    
    phi = lat + psi .* cos(az*pi);

    phi(phi > 0.416)  =  0.416;
    phi(phi < -0.416) = -0.416;
    
    % geodetic longitude of the earth projection of the ionospheric intersection point
    lambda = lon + ((psi.*sin(az*pi)) ./ cos(phi*pi));

    % geomagnetic latitude of the earth projection of the ionospheric intersection point
    ro = phi + 0.064*cos((lambda-1.617)*pi);

    % local time in seconds
    t = lambda*43200 + time;

    if t >= 86400
        t = t - 86400;
    elseif t < 0
        t = t + 86400;
    end
    
    % Obliquity factor
    f = 1 + 16*(0.53-el).^3;     

    a = a0 + a1*ro + a2*ro.^2 + a3*ro.^3;
    a(a < 0) = 0;
    
    p = b0 + b1*ro + b2*ro.^2 + b3*ro.^3;
    p(p < 72000) = 72000;
    
    x = (2*pi*(t-50400)) ./ p;

    %ionospheric delay
    index = find(abs(x) < 1.57);
    delay(index,1) = c * f(index) .* (5e-9 + a(index) .* (1 - (x(index).^2)/2 + (x(index).^4)/24));
    
    index = find(abs(x) >= 1.57);
    delay(index,1) = c * f(index) .* 5e-9;
end
