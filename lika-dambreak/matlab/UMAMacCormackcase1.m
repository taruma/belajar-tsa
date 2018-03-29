% Program Saint-Venant
% MacCormack Scheme
% Case: DAM BREAK
% Script by Taruma (25017046), hi@taruma.info
% Documentation: Hidrolika T04 20180222.pdf
% Github this project: https://github.com/taruma/belajar-tsa
% Folder: lika-dambreak/matlab

% Initiating Program
clc;    % Clear command window
clf;    % Clear current figure window
clear;  % Remove items from workspace, freeing up system memory

% Note:
% x Axis define by i
% time define by t
% Assume channel cross-section is trapesium

% Defining Constant
% Start---
g = 9.8;    % Gravity Acceleration
n = 0;      % Manning Slope
beta = 1;   % Correction

dx = 0.5;   % Grid length
dt = 0.1;   % Time step

% Secant Method 
dy = 0.01;  % Delta y for secant solution
yit = 10;   % Iteration

time = 50; % Duration simulation in seconds
tout = 1;   % Plotting each 'tout' second(s)

% Water Elevation
h1 = 1.5;   % Upstream head in (meter)
h2 = 1;     % Downstream head in (meter)

% Channel Properties
lf = 55;    % Length of channel in (meter)
m = 0;      % Slope of trapesium cross-section

% Grid/Time Step
imax = lf/dx+1;   % Total grid
tmax = time/dt; % Time Step (iteration)

% Width of channel (define for each grid)
set_b = 1;
for i = 1:imax
    b(i) = set_b;
end
% ---End

% Defining Variable and Initial Condition
% Start---
% Assigning variable as 0 value in n-dimension array
Q = zeros(1, imax); % Discharge Q = Q0_0
z = zeros(1, imax); % Elevation of z, Assigning 0
A = zeros(1, imax); % Channel Area
h = zeros(1, imax); % Depth
R = zeros(1, imax); % Hydraulic Radius
P = zeros(1, imax); % Wet Perimeter

% Predictor Variable (Only for predictor step), using p-endfix
Qp = Q;
zp = z;
Ap = A;
hp = h;

% For analytical solution
hex = h;

% Initial Condition
for i = 1:imax
    if i <= imax/2
        h(i) = h1;
    else
        h(i) = h2;
    end
    A(i) = (b(i) + m * h(i)) * h(i);
    P(i) = b(i) + 2*h(i)*sqrt(1 + m.^2);
    R(i) = A(i) / P(i);
end

% Variable depend by time
An = A; % An = Ap1_0
Qn = Q; % Qn = Qp1_0
% ---End

% Start PROGRAM
% Iteration by time
p3 = plot(z, 'black');
for t = 1:tmax
    % Start of Predictor Step----
    for i = 2:imax-1
        % Continuity-Initial
        A(i) = (b(i) + m*h(i))*h(i);
        P(i) = b(i) + 2*h(i)*sqrt(1+m^2);
        R(i) = A(i)/P(i);
        % Continuity-Predictor
        Ap(i) = A(i) - dt/dx * (Q(i+1)-Q(i));
        
        % Momentum-Initial
        I2p = beta * ((Q(i+1)^2/A(i+1)) - (Q(i)^2/A(i)))/dx;
        I3p = g * A(i) * ((h(i+1)+z(i+1))-(h(i)+z(i)))/dx;
        I4p = g * Q(i) * abs(Q(i)) * n^2 / (A(i)*R(i)^(4/3));
        % Momentum-Predictor
        Qp(i) = Q(i) - dt * (I2p + I3p + I4p);        
    end
    
    % Boundary Condition-Predictor
    Qp(1) = Qp(2);
    Qp(imax) = Qp(imax-1);
    Ap(1) = Ap(2);
    Ap(imax) = Ap(imax-1);
    
    % Find depth (h) value
    for i = 1:imax
        % Solve for y by Secant Method
        yawal = 5;
        for j = 1:yit 
            yplus = yawal + dy;
            ymin = yawal - dy;
            
            Fawal = Ap(i) - b*yawal - m * yawal^2;
            Fplus = Ap(i) - b*yplus - m * yplus^2;
            Fmin = Ap(i) - b*ymin - m * ymin^2;
            
            dF = (Fplus - Fmin) / (2*dy);

            yfinal = yawal - Fawal/dF;
            yawal = yfinal;
        end
        hp(i) = yfinal;
        zp(i) = z(i);
    end
    % ---End of Predictor Step
    
    % Start of Corrector Step ---
    for i = 2:imax-1
        % Continuity-Initial
        A(i) = (b(i) + m*h(i))*h(i);
        P(i) = b(i) + 2*h(i)*sqrt(1+m^2);
        R(i) = A(i)/P(i);        
        % Continuity-Corrector
        Aph = (A(i) + Ap(i))/2;
        An(i) = Aph - dt/(2*dx) * (Qp(i)-Qp(i-1));        
        
        % Momentum-Initial
        I2c = beta * ((Qp(i)^2/Ap(i)) - (Qp(i-1)^2/Ap(i-1)))/dx;
        I3c = g * A(i) * ((hp(i)+zp(i))-(hp(i-1)+zp(i-1)))/dx;
        I4c = g * Q(i) * abs(Q(i)) * n^2 / (A(i)*R(i)^(4/3));  
        % Momentum-Corrector
        Qph = (Q(i)+Qp(i))/2;
        Qn(i) = Qph - dt/2 * (I2c + I3c + I4c);        
    end
    
    % Boundary Condition-Corrector/Final
    Qn(1) = Qn(2);
    Qn(imax) = Qn(imax-1);
    An(1) = An(2);
    An(imax) = An(imax-1);    
    % ---End of Corrector Step
    
    % Current -> Next Step
    A = An;
    Q = Qn;
    
    % Find depth with Current
    for i = 1:imax
        yawal = 5;
        for j = 1:yit % 5 iterasi
            yplus = yawal + dy;
            ymin = yawal - dy;
            
            Fawal = A(i) - b*yawal - m * yawal^2;
            Fplus = A(i) - b*yplus - m * yplus^2;
            Fmin = A(i) - b*ymin - m * ymin^2;
            
            dF = (Fplus - Fmin) / (2*dy);

            yfinal = yawal - Fawal/dF;
            yawal = yfinal;
        end       
        h(i) = yfinal;
    end    
    
    % --- End of Numerical Solution
    
    % Analytic Solution ---
    hm = 1.237805;
    cm = sqrt(g * hm);
    xa = (imax - 1)/2 * dx - (t * dt) * sqrt(g*1.5);
    xb = (imax - 1)/2 * dx + (t * dt) * (2 * sqrt(g * 1.5) - 3 * cm);
    xc = (imax - 1)/2 * dx + (t * dt) * (2 * cm^2 * (sqrt(g * 1.5) - cm) / (cm^2-g*1));
    
    for i = 1 : imax
        if (i - 1) * dx <= xa
            hex(i) = 1.5;
        elseif (i-1) * dx <= xb
            hex(i) = 4 / (9*g)*(sqrt(g*1.5)-(i*dx-imax/2*dx)/(2*(t*dt)))^2;
        elseif and((i-1)*dx > xb, (i-1) * dx <= xc)
            hex(i) = hm;
        else
            hex(i) = 1;
        end
    end
    
    % Graph Generator
    if mod(t, tout) == 0
        p1 = plot(h,'red');
        hold on

        p2 = plot(hex,'blue');
        hold on
        
        axis manual
        axis ([1 imax 0 2]);
        title('Dam Break MacCormack Method');
        xlabel('Grid (x)');
        ylabel('Elevasi(m)');
        pause(0.001)
        hold off
    end
end