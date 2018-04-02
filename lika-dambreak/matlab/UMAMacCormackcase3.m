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

% ---Defining Constant and something important
g = 9.81;   % Gravity Acceleration
n = 0;      % Manning Slope
beta = 1;   % Correction
dx = 0.5;   % Grid length
dt = 0.1;   % Time step
            % --Secant Method 
dy = 0.1;   % Delta y for secant solution
yit = 10;   % Iteration
            % --Time
time = 60;  % Duration simulation in seconds
tone = 1/dt;        % (tone) iteration for one second
tout = tone*0.01;   % Plotting each 'tout' second(s)
            % --Water Elevation
h1 = 1.5;   % Upstream head in (meter) elevation
h2 = 1;     % Downstream head in (meter) elevation
            % --Channel Properties
lf = 55;    % Length of channel in (meter)
lc = 43;    % Length of contraction in (meter) from left
sc = lc;    % Slope change at (lc) meter from left
m = 1;      % m of trapesium cross-section

                % Grid/Time Step
imax = lf/dx+1; % Total grid
tmax = time/dt; % Time Step (iteration)

                    % Boundary Condition Variable
wall_check = false; % True for Wall Condition

                    % ---Assigning width value
b1 = 10;            % Width b1
b2 = b1;            % Width b2
b = zeros(1, imax); % Preset b array
ic = lc/dx+1;       % Start at (ic), width is change
for i = 1:imax
    if i <= ic
        b(i) = b1;
    else
        b(i) = b2;
    end
end
                    % Defining Variable and Initial Condition
                    % Assigning variable as 0 value in n-dimension array
Q = zeros(1, imax); % Discharge Q = Q0_0
z = zeros(1, imax); % Elevation of z, Assigning 0
A = zeros(1, imax); % Channel Area
h = zeros(1, imax); % Depth
R = zeros(1, imax); % Hydraulic Radius
P = zeros(1, imax); % Wet Perimeter

% Reasign p value (elevation of base channel) with to 0 x-axis
chan_ev = zeros(1, imax);
chan_up = 0;
slope1 = 0;         % '-' mean ascending
slope2 = slope1;         % 
chan_ev(1) = chan_up;   % start elevation
is = sc/dx+1;           % 
for i = 2:imax
    if i <= is
        chan_ev(i) = chan_ev(i-1) - slope1*dx;
    else
        chan_ev(i) = chan_ev(i-1) - slope2*dx;
    end
end
% Reasign z value, shift datum to never negative position (lowest)
lowest = min(chan_ev);          % finding the lowest value
for i = 1:imax
    z(i) = chan_ev(i)-lowest;   % changing datum (never negative)
end
            % Another Something.. 
Qp = Q;     % Predictor Variable (Only for predictor step), using p-endfix
zp = z;
Ap = A;
hp = h;
Rp = R;
Pp = P;
hex = h;    % For analytical solution
            % Position in meter
pos_x = zeros(1, imax);
for i = 1:imax
    pos_x(i) = i*dx;
end
val_x = 1:imax;
% ---End something

% Initial Condition
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

% Checking Value (temporary)
% h
% z
% h+z
% z(1)
% z(imax)

% Variable depend by time
An = A; % An = Ap1_0
Qn = Q; % Qn = Qp1_0
% ---End

% Start PROGRAM
% Iteration by time
highest = 0;
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
    if wall_check
        Qp(1) = 0;
        Qp(imax) = 0;
    else
        Qp(1) = Qp(2); % Ignore the warning message
        Qp(imax) = Qp(imax-1);
    end
    Ap(1) = Ap(2);
    Ap(imax) = Ap(imax-1);
    
    % Find depth (h) value
    for i = 1:imax
        % Solve for y by Secant Method
        yawal = 5;
        for j = 1:yit 
            yplus = yawal + dy;
            ymin = yawal - dy;
            
            Fawal = Ap(i) - b(i)*yawal - m * yawal^2;
            Fplus = Ap(i) - b(i)*yplus - m * yplus^2;
            Fmin = Ap(i) - b(i)*ymin - m * ymin^2;
            
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
        Ap(i) = (b(i) + m*hp(i))*hp(i);
        Pp(i) = b(i) + 2*hp(i)*sqrt(1+m^2);
        Rp(i) = Ap(i)/Pp(i);        
        % Continuity-Corrector
        Aph = (A(i) + Ap(i))/2;
        An(i) = Aph - dt/(2*dx) * (Qp(i)-Qp(i-1));        
        
        % Momentum-Initial
        I2c = beta * ((Qp(i)^2/Ap(i)) - (Qp(i-1)^2/Ap(i-1)))/dx;
        I3c = g * Ap(i) * ((hp(i)+zp(i))-(hp(i-1)+zp(i-1)))/dx;
        I4c = g * Qp(i) * abs(Qp(i)) * n^2 / (Ap(i)*Rp(i)^(4/3));  
        % Momentum-Corrector
        Qph = (Q(i)+Qp(i))/2;
        Qn(i) = Qph - dt/2 * (I2c + I3c + I4c);        
    end
    
    % Boundary Condition-Corrector/Final
    if wall_check
        Qn(1) = 0;
        Qn(imax) = 0;
    else
        Qn(1) = Qn(2); % Ignore the warning message
        Qn(imax) = Qn(imax-1);
    end
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
            
            Fawal = A(i) - b(i)*yawal - m * yawal^2;
            Fplus = A(i) - b(i)*yplus - m * yplus^2;
            Fmin = A(i) - b(i)*ymin - m * ymin^2;
            
            dF = (Fplus - Fmin) / (2*dy);

            yfinal = yawal - Fawal/dF;
            yawal = yfinal;
        end       
        h(i) = yfinal;
    end    
    
    % --- End of Numerical Solution
    
    % Assigning lowest/highest value to graph
    newlow = min(h);
    newhigh = max(h);
    
    if newlow < lowest
        lowest = newlow;
    elseif newhigh > highest
        highest = newhigh;
    end
    
    
    % Analytic Solution ---
    hm = 1.23;
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
        p1 = plot(val_x,h + chan_ev,'red');
        hold on
        p3 = plot(val_x,chan_ev, 'black');
        hold on
        p2 = plot(hex,'blue');
        hold on
        
        %xticks(pos_x);

        axis manual
        axis ([1 imax lowest-0.1 highest+0.1]);
        title('Dam Break MacCormack Method');
        xlabel('Grid');
        ylabel('Elevasi(m)');
        pause(0.001)
        hold off
    end
end