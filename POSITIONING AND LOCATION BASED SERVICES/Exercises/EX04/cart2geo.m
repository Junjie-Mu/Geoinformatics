function [lat, lon, h] = cart2geo(X,Y,Z)
        a = 6378137;
        b = 6356752.314;
        eb = sqrt((a^2-b^2)/b^2);
        r = sqrt(X^2+Y^2);
        e = sqrt(1 - b^2/a^2);
   
        Psi = atan2((Z*a),(r*b));
        Lambda = atan2(Y,X);
        Fai = atan2((Z+eb^2*b*sin(Psi)^3),(r-e^2*a*cos(Psi)^3));
        Rn = a / sqrt(1-e^2*sin(Fai)^2);
        %black magic!
        h = r/cos(Fai);
        lat = Fai/pi*180;
        lon = Lambda/pi*180;
end