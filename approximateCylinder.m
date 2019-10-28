function [points] = approximateCylinder(n, r)
    % Ladies & Gentlemen, Will Crain
    % How does he do it!?
    
    theta = 0:2*pi/n:2*pi;
    thetaOffset = pi/n;
    
    points = fliplr([r*cos(theta(1:end-1) + thetaOffset); r*sin(theta(1:end-1) + thetaOffset)]');
end