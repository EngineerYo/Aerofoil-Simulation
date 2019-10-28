function [V] = getVelocity(sPanel, tPanel)
    deltR = sPanel(1, 13:14) - tPanel(1, 13:14);
    RdotT = dot(deltR, tPanel(1, 11:12));
    Lb = sPanel(1, 7);
    dR2rDt2 = sqrt(norm(deltR)^2 - RdotT^2);
    
    V = ((deltR -  tPanel(1, 11:12)*RdotT)/(2*pi*dR2rDt2)) * (satan((Lb/2 - RdotT), dR2rDt2) - satan((-Lb/2 - RdotT), dR2rDt2)) + ...
        -(tPanel(1, 11:12)/(4*pi)) * (log(1+((Lb/2 - RdotT)/dR2rDt2)^2) - log(1+((-Lb/2 - RdotT)/dR2rDt2)^2));
end