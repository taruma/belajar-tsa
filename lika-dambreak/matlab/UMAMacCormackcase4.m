%% INFORMATION! READ ME FIRST
% Program Saint-Venant: MacCormack Scheme
% Case: 4 - Tidal Wave 
% Script by Taruma (25017046), hi@taruma.info
% Documentation: N/A
% Github this project: https://github.com/taruma/belajar-tsa
% Folder: lika-dambreak/matlab
% WARNING, this script is using Function in script, Use 2016b or later!!

%% --- Initiating Program
close all;  % Close all window
clc;        % Clear command window
clear;      % Remove items from workspace, freeing up system memory

%% --- Define Variable and Constant
g = 9.81;   % Gravity Acceleration
n = 0;      % Manning Slope
beta = 1;   % Correction

%% --- Channel (Elevation and Properties)
% --Water Elevation
h1 = 2;     % Upstream head in (meter) elevation
h2 = h1;    % Downstream head in (meter) elevation

% --Channel Properties
lf = 100;   % Length of channel in (meter)
m = 1;      % m of trapesium cross-section

%% --- Secant Method
yawal = 5;  % Guess value for y
dy = 0.1;   % Delta y
yit = 10;   % Iteration

%% --- Grid and Time
dx = 0.5;   % Grid length
dt = 0.1;   % Time step
% -- Time
time = 100;         % Duration simulation in seconds
tone = 1/dt;        % (tone) iteration for one second
tout = tone*0.1;    % Plotting each 'tout' second(s)
% Grid/Time Step
imax = lf/dx+1; % Total grid
tmax = time/dt; % Time Step (iteration)

%% --- Tidal properties
periode = 5;
ampli = 0.1;

%% --- Width and Position
% INPUT
lc = 80;    % Length of contraction in (meter) from left
b1 = 10;     % Width b1
b2 = b1;    % Width b2

% -- Assigning width value
b = zeros(1, imax); % Preset b array
ic = lc/dx+1;       % Start at (ic), width is changing
for i = 1:imax
    if i <= ic
        b(i) = b1;
    else
        b(i) = b2;
    end
end

% - Assigning position value & absis value
pos_x = zeros(1, imax);
for i = 1:imax
    pos_x(i) = i*dx;
end
val_x = 1:imax;         % Absis value

%% --- Channel and Datum Elevation
% Input
chan_up = -0.3;               % Channel elevation at grid 1 (left)
slope1 = -0.3/lf;                 % '-' mean ascending
slope2 = 0;                 % 
sc = 0.9*lf;    % Slope change at (lc) meter from left

% Channel Elevation
chan_ev = zeros(1, imax);   % Preset chan_ev
chan_ev(1) = chan_up;       % start elevation
is = sc/dx+1;               % at grid (is), slope is changing
for i = 2:imax
    if i <= is
        chan_ev(i) = chan_ev(i-1) - slope1*dx;
    else
        chan_ev(i) = chan_ev(i-1) - slope2*dx;
    end
end

% Reasign z value, shift datum to never negative position (lowest)
z = zeros(1, imax);             % Elevation of z, Assigning 0
lowest = min(chan_ev);          % finding the lowest value
for i = 1:imax
    z(i) = chan_ev(i)-lowest;   % changing datum (never negative)
end

%% Misc.. (wall_check, ...)
% Boundary Condition Variable
wall_check = true;  % True for Wall Condition
highest = 0;        % Highest value (for graph boundary)

%% Defining Variable and Initial Condition
% This code for preset dimension of variable
% Assigning variable as 0 value in n-dimension array
Q = zeros(1, imax); % Discharge Q = Q0_0
A = zeros(1, imax); % Channel Area
h = zeros(1, imax); % Depth
R = zeros(1, imax); % Hydraulic Radius
P = zeros(1, imax); % Wet Perimeter

% Define Predictor
Qp = Q;     % Predictor Variable (Only for predictor step), using -p prefix
zp = z;
Ap = A;
hp = h;

% Variable depend by time
An = A; % An = Ap1_0
Qn = Q; % Qn = Qp1_0

%% PROG: INITIAL CONDITION
% Assign water elevation and Calc. channel prop. 
for i = 1:imax
    if i <= imax/2
        h(i) = h1 - chan_ev(i);
    else
        h(i) = h2 - chan_ev(i);
    end
    A(i) = (b(i) + m * h(i)) * h(i);
    P(i) = b(i) + 2*h(i)*sqrt(1 + m.^2);
    R(i) = A(i) / P(i);
end
% ---End

%% PROG: MACCORMACK SCHEME
% Iteration by time
for t = 1:tmax
    tnow = t*dt;                    % Calc. (real) time
    
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
    if wall_check % IGNORE THE WARNING MESSAGE
        Qp(1) = 0; 
    else
        Qp(1) = Qp(2);
    end
    Qp(imax) = Qp(imax-1); 
    Ap(1) = Ap(2);
    hp(imax) = h2 + tidalwave(tnow, periode, ampli);
    hp(imax) = hp(imax) - chan_ev(imax); 
    Ap(imax) = (b(imax) + m*hp(imax))*hp(imax);

    zp = z;
    % Find depth (h) value
    for i = 1:imax-1
        hp(i) = secantmet(yawal, dy, yit, [Ap(i) b(i) m]);
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
    if wall_check % IGNORE THE WARNING MESSAGE
        Qn(1) = 0; 
    else
        Qn(1) = Qn(2);
    end
    
    Qn(imax) = Qn(imax-1);
    An(1) = An(2);
    h(imax) = h2 + tidalwave(tnow, periode, ampli); 
    h(imax) = h(imax) - chan_ev(imax);
    An(imax) = (b(imax) + m*h(imax))*h(imax);
    % ---End of Corrector Step
    
    % Current -> Next Step
    A = An;
    Q = Qn;
    
    % Find depth with Current
    for i = 1:imax-1
        h(i) = secantmet(yawal, dy, yit, [A(i) b(i) m]);
    end    
    
    % --- End of Numerical Solution
    
    % Assigning lowest/highest value to graph
    newlow = min(h);
    newhigh = max(h);
    
    % Check Error and stop the program if conditional are met
    if newhigh >= 100 || newlow <= -100
        return
    end
    
    if newlow < lowest
        lowest = newlow;
    elseif newhigh > highest
        highest = newhigh;
    end
    
    % Graph Generator
    if mod(t, tout) == 0
        subplot(2,1,1);
        p1 = plot(val_x,h + chan_ev,'red');
        hold on
        p3 = plot(val_x,chan_ev, 'black');
        hold on
        axis manual
        axis ([1 imax lowest-0.1 highest+0.1]);
        title('Tidal Wave MacCormack Method');
        %xlabel('Grid');
        ylabel('Elevasi(m)');
        pause(0.001)
        hold off

        subplot(2,1,2);
        p2 = plot(val_x,Q,'green');
        hold on
        title('Q');
        axis manual
        axis ([1 imax -30 30]);
        xlabel('Grid');
        ylabel('Q (m3/s)');
        pause(0.001)
        hold off

    end
end

%% FUNCTION: PROGRAM
function h = tidalwave(t, T, A)
    % t = Time
    % T = Period
    % A = Amplitude
    h = sin(2*pi*t/T)*A;
end

function y = secantmet(yawal, dy, yit, Abm)
%        Example Code:
%         yawal = 5;
%         for j = 1:yit % 5 iterasi
%             yplus = yawal + dy;
%             ymin = yawal - dy;
%             
%             Fawal = A(i) - b(i)*yawal - m * yawal^2;
%             Fplus = A(i) - b(i)*yplus - m * yplus^2;
%             Fmin = A(i) - b(i)*ymin - m * ymin^2;
%             
%             dF = (Fplus - Fmin) / (2*dy);
% 
%             yfinal = yawal - Fawal/dF;
%             yawal = yfinal;
%         end       
%        h(i) = yfinal;
    A = Abm(1);
    b = Abm(2);
    m = Abm(3);
    for i = 1:yit
        yplus = yawal + dy;
        ymin = yawal - dy;
        
        Fawal = A - b*yawal - m * yawal^2;
        Fplus = A - b*yplus - m * yplus^2;
        Fmin =  A - b*ymin  - m * ymin^2;

        dF = (Fplus - Fmin) / (2*dy);

        yfinal = yawal - Fawal/dF;
        yawal = yfinal;
    end
    y = yfinal; 
end