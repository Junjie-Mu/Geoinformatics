function [tropoDelay] = tropo_correction(h, el)

%--------------------------------------------------------------------------
%  Saastamoinen model
%  Computation of the pseudorange correction due to tropospheric refraction.
%    
%   Parameters
%   ----------
%   h : height f the receiver [m].
%   el : elevation of the satellite with respect the receiver [degrees].
% 
%   Returns
%   -------
%   tropoDelay : tropospheric effect [m].
% 
%--------------------------------------------------------------------------

if (h < 5000)
    %conversion to radians
    el = abs(el) * pi/180;
    
    %Standard atmosphere - Berg, 1948 (Bernese)
    Po = 1013.25;                                 % pressure [mbar]
    To = 291.15;                                  % temperature [K]
    Ho = 50/100;                                  % humidity [%]
    ho = 0;
    
    height = h - ho;                              % h is the ellipsoidal height of the receiver
    
    Pr = Po * (1-0.0000226 * height)^5.225; 
    Tr = To - 0.0065 * height;
    Hr = Ho * exp(-0.0006396 * height);
    
    er = Hr .* exp(-37.2465 + 0.213166*Tr - 0.000256908*Tr.^2);
    tropoDelay = (0.002277 ./ sin(el)) .* (Pr + (1255./Tr + 0.05).*er - (tan(el)).^-2);
else
    tropoDelay = 0;
end
end

%Marianna Alghisi