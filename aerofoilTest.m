%% Crain, William
% Aerofoil example
% September 25, 2018
% ME 482

close all; clc; clear all; format compact;

%% Constants

% R = importdata('airfoil2412_32.txt');
R = [-sqrt(3), 0; 0, 1; 0, -1];

%% Work

cells = zeros(size(R, 1), 16);
% Cells
    % 1         Start node
    % 2         End node
    % 3, 4      Start X, Start Y
    % 5, 6      End X, End Y
    % 7         Length
    % 8         Angle with X axis
    % 9, 10     Normal X, Normal Y
    % 11, 12    Tangent X, Tangent Y
    % 13, 14    Midpoint X, Midpoint Y
    % 15, 16    u�V, v�V

for i=1:size(R, 1)
    target = i+1;
    if i==size(R, 1)
        target = 1;
    end
    
    dx = R(target, 1) - R(i, 1);
    dy = R(target, 2) - R(i, 2);
    
    length = sqrt(dx^2 + dy^2);
    
    angle = atan(dy/dx);
    if dy/dx < 0
        angle = angle + pi;
    end
    
    nX = -sin(angle);
    nY = cos(angle);
    
    tX = cos(angle);
    tY = sin(angle);
    
    mX = R(i, 1) + dx/2;
    mY = R(i, 2) + dy/2;
    
    cells(i, :) = round([i, target, R(i, 1), R(i, 2), R(target, 1), R(target, 2), length, angle, nX, nY, tX, tY, mX, mY, 0, 0], 5);
end

for i = 1:size(R, 1)
    
    Vi = 0;
    
    for v = 1:size(R, 1)
        if (v ~= i)
            sPanel = cells(i, :);
            tPanel = cells(v, :);
            
            deltRi = tPanel(1, 13:14) - sPanel(1, 13:14);
            RdotT = dot(sPanel(1, 13:14), sPanel(1, 11:12));
            Lb = sPanel(1, 7);
            deltR2 = norm(tPanel(1, 13:14))^2;
            dR2rDt2 = sqrt(deltR2 - RdotT^2);
            
            p1 = (deltRi - tPanel(1, 11:12)*RdotT)/(2*pi*dR2rDt2^2);
            p2 = atan((Lb/2 - RdotT)/dR2rDt2) - atan((-Lb/2 - RdotT)/dR2rDt2);
            p3 = -tPanel(1, 11:12)/(4*pi);
            p4 = log(1+((Lb/2 - RdotT)/dR2rDt2)^2) - log(1+((-Lb/2 - RdotT)/dR2rDt2)^2);
            
            Vi = Vi + p1*p2 + p3*p4;            
        end
        cells(i, 15:16) = Vi;
    end
end

tab = array2table(cells, 'VariableNames', {'P0', 'P1', 'P0x', 'P0y', 'P1x', 'P1y', 'L', 'Angle', 'Nx', 'Ny', 'Tx', 'Ty', 'Mx', 'My', 'Vx', 'Vy'});

tab(:, 15:16)

%%

hold on
plot(cells(:, 3:2:5), cells(:, 4:2:6), 'k')
% axis([-0.5, 1.5, -1, 1])
% axis([-2, 1, -1.5, 1.5])