%% Crain, William
% Template
% November 5, 2018
% ME 482

clc; clear variables; format compact;
% close all;

%% Aerofoils

% Test Aerofoil
% R = importdata('airfoil2412_32.txt');

% Test Aerofoil 2
% R = importdata('airfoil6409.txt');
% R = importdata('NACA4415.txt');

% Triangle aerofoil
R = [-sqrt(3), 0; 0, 1; 0, -1];

% Cylinder aerofoil
% R = approximateCylinder(32, 1);

% Flat aerofoil
% R = [-1, 0; 1, 0];

R = transformPoints(R);

%% Initial setup

cells = zeros(size(R, 1), 16);
alpha = 0;
U = 1;
gamma = 1;

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
    % 15, 16    B and Lambda
    % 17        Pressure coefficient

for i=1:size(R, 1)
    
    % Account for i = size(R, 1), target 1st vertex for enclosed shape
    target = i+1;
    if i==size(R, 1)
        target = 1;
    end
    
    dx = R(target, 1) - R(i, 1);
    dy = R(target, 2) - R(i, 2);
    
    length = sqrt(dx^2 + dy^2);
    
    angle = atan2(dy, dx);
    
    nX = -sin(angle);
    nY = cos(angle);
    
    tX = cos(angle);
    tY = sin(angle);
    
    mX = R(i, 1) + dx/2;
    mY = R(i, 2) + dy/2;
    
    B = -U * sin(alpha - angle);
    
    cells(i, :) = [i, target, R(i, 1), R(i, 2), R(target, 1), R(target, 2), length, angle, nX, nY, tX, tY, mX, mY, B, 0];
end

% targetBound accounts for flat plates
% For #Points >= 3, they're #Points panels
% For #Points < 3, they're #Points-1 panels
targetBound = size(R, 1);
if (size(R, 1) < 3)
    targetBound = size(R, 1) - 1;
end

%% Calculate Velocity

V = zeros(targetBound, targetBound, 4);

for i = 1:targetBound    
    for v = 1:targetBound
        
        sPanel = cells(i, :);
        tPanel = cells(v, :);
        
        % If we target panel B while calculating panel B,
        % Make target panel equal to panel B, but slightly
        % offset in the counter clockwise direction
        if (i == v)
            tPanel = sPanel;
            
            % Add a small number such that the normal and 
            % tangent vectors point in the correct directions
            dx = -sPanel(1, 9)/1e9;
            dy = -sPanel(1, 10)/1e9;
            
            tPanel(1, 13) = tPanel(1, 13) + dx;
            tPanel(1, 14) = tPanel(1, 14) + dy;
        end
        
        % getVelocity is a function that calculates velocity
        V(i, v, 1:2) = getVelocity(sPanel, tPanel);
        V(i, v, 3:4) = [V(i, v, 2), -V(i, v, 1)];
    end
end

%% Calculate sources

n = zeros(size(V, 1), size(V, 1), 2);
nv = zeros(size(V, 1), size(V, 1), 1);

for i = 1:targetBound
    for v = 1:targetBound
        n(i, v, :) = cells(i, 9:10);
        nv(i, v) = dot(V(i, v, 1:2), n(i, v, :));
    end
end

% This calculates lambda
cells(:, 16) = cells(:, 15)' / nv;

%% V = Vab * lambda + Uinfinity

numCells = 20;
offset = 5;

xmin = min(R(:, 1))-offset;
xmax = max(R(:, 1))+offset;
xstep = (xmax-xmin)/numCells;
x = xmin:xstep:xmax;

ymin = min(R(:, 2))-offset;
ymax = max(R(:, 2))+offset;     
ystep = (ymax-ymin)/numCells;   
y = ymin:ystep:ymax;

Rf = zeros(numCells, numCells, 2);
for i = 1:1:numCells
    for v = 1:1:numCells
        Rf(i, v, :) = [x(i), y(v)];
    end
end

Vplot = zeros(numCells, numCells, 2);
for i = 1:1:numCells
    for v = 1:1:numCells
        for t = 1:1:targetBound
            tPanel = cells(t, :);
            sPanel = tPanel;
            sPanel(1, 13) = x(i);
            sPanel(1, 14) = y(v);
            
            vB = getVelocity(sPanel, tPanel) * tPanel(1, 16);
            
            Vplot(i, v, 1) = Vplot(i, v, 1) + vB(1);
            Vplot(i, v, 2) = Vplot(i, v, 2) + vB(2);
        end
        Vplot(i, v, 1) = Vplot(i, v, 1) + U*cos(alpha) - y(v)*gamma/sqrt(x(i)^2+y(v)^2);
        Vplot(i, v, 2) = Vplot(i, v, 2) + U*sin(alpha) + x(i)*gamma/sqrt(x(i)^2+y(v)^2);
    end
end

%% Find pressure coefficients and tangential velocity

Vt = zeros(1, size(R, 1));

for i = 1:size(V, 1)
    sPanel = cells(i, :);
    
    Vt(1, i) = U*cos(alpha-sPanel(1, 8)) + dot(sPanel(1, 11:12), [V(i, i, 1), V(i, i, 2)]);
    for v = 1:size(V, 2)
        tPanel = cells(v, :);
        Vt(1, i) = Vt(1, i) + dot(sPanel(1, 11:12), [V(i, v, 1), V(i, v, 2)])*tPanel(1, 16) + gamma/(2*pi*sqrt(sPanel(1, 13)^2 + sPanel(1, 14)^2));
    end
end

cells(:, 17) = 1 - (Vt(1, :)/U).^2;

%% Clean up velocities

[X, Y] = meshgrid(x, y);
[in, on] = inpolygon(X, Y, R(:, 1), R(:, 2));

for i = 1:numCells
    for v = 1:numCells
        if (in(v, i) == 1 || on(v, i) == 1)
            Vplot(i, v, :) = [0, 0];
        end
    end
end

%% Plotting

figure
hold on
title('Gamma of 0.8')
plot(cells(:, 3:2:5), cells(:, 4:2:6), 'k', 'LineWidth', 2);
axis([xmin, xmax, ymin, ymax])
quiver(Rf(:, :, 1), Rf(:, :, 2), Vplot(:, :, 1), Vplot(:, :, 2), 'r', 'LineWidth', 1);
pbaspect([1, 1, 1]);

figure
hold on
scatter(1:size(R, 1), cells(:, 17), 'kx')
plot(1:size(R, 1), cells(:, 17), 'k', 'LineWidth', 0.5)
title('Pressure coefficients')
ticks = 0:(size(R, 1)/10):size(R, 1);
xticks(ticks)
xticklabels(ticks*360/size(R, 1))
pbaspect([1, 1, 1]);

figure
hold on
scatter(1:size(R, 1), Vt(1, :), 'kx')
plot(1:size(R, 1), Vt(1, :), 'k', 'LineWidth', 0.5)
title('Tangential Velocity')
ticks = 0:(size(R, 1)/10):size(R, 1);
xticks(ticks)
xticklabels(ticks*360/size(R, 1))
pbaspect([1, 1, 1]);

%% To-do list

% =)

%% Cleanup

% Turn 'cells' into a table, which has labeled collumns
cells = array2table(cells, 'VariableNames', {'P0', 'P1', 'P0x', 'P0y', 'P1x', 'P1y', 'L', 'Angle', 'Nx', 'Ny', 'Tx', 'Ty', 'Mx', 'My', 'B', 'Lambda', 'Cp'});

% Clean up memory by clearing out intermediate calculated values
toClean = {'angle', 'gamma', 'deltR', 'dR2rDt2', 'dx', 'dy', 'i', 'Lb', 'length', 'mX', 'mY', 'nX', 'nY', 'tX', 'tY', 'p1', 'p2', 'p3', 'p4', 'RdotT', 'sPanel', 'target', 'targetBound', 'tPanel', 'v', 'Vi', 'toClean', 'lambda', 'sign', 'alpha', 'ans', 'ldX', 'ldY', 'size', 't', 'U', 'vt', 'vz', 'xmax', 'xmin', 'xstep', 'ymax', 'ymin', 'ystep', 'vB', 'Rf', 'Rf', 'numCells', 'A', 'B', 'C', 'D', 'E', 'I', 'J', 'n', 'offset', 'Sj', 'VdotT', 'Vt', 'X', 'Y', 'x', 'y', 'on', 'in', 'ticks'};
clear(toClean{:})