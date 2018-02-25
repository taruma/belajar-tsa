% Program Saint-venant
% Kasus: Dam Break
% Leap Frog Scheme and Lax Wendroff
% Script by Taruma (25017046), hi@taruma.info


% Initiate Program
clc;
clear;

% Sumbu x didefinisikan dengan i
% Waktu didefinisikan dengan t
% Saluran diasumsikan berbentuk segi-empat

% ===============================
% Peubah dan Konstanta 
g = 9.8;    % percepatan gravitasi
n = 0;      % koefisien kekasaran
beta = 1;   % koefisien koreksi

dx = 0.5;   % length/grid step
dt = 0.0005;  % time step

time = 7;   % Waktu simulasi
tout = 100;  % Plot tiap waktu tout

% Kondisi Dam
b = 1;      % Lebar saluran
h1 = 1.5;   % Tinggi muka air di kiri (hulu)
h2 = 1;     % Tinggi muka air di kanan (hilir)
lf = 54;    % Panjang saluran

% Grid
imax = lf/dx;   % Jumlah grid untuk jarak (panjang saluran/jarak grid)
tmax = time/dt; % Jumlah time step (waktu simulasi/time step)

% Pendefinisian peubah dan kondisi inisial
% Peubah nol
% ------ Leap Frog --------
Q = zeros(1, imax); % Debit Q(t,i) = Q; Q(t+1,i) = Qn; Q(t-1,i) = Qo
z = zeros(1, imax); % tinggi datum terhadap dasar permukaan saluran
A = zeros(1, imax); % Luas permukaan
R = zeros(1, imax); % Jari-jari hidrolis = Luas permukaan / Keliling basah
h = zeros(1, imax); % kedalaman air
P = zeros(1, imax); % Keliling basah

% ------ Solusi Analitik -----
hex = h;            % kedalaman air exact (solusi analitik)

% ------ Lax-Wendroff --------
Qlw = Q;
zlw = z;
Alw = A;
Rlw = R;
hlw = h;            % kedalaman air lax-wendroff
Plw = P; 

%Initial Condition
for i = 1:imax
    if i <= imax/2
        h(i) = h1;
    else
        h(i) = h2;
    end
    A(i) = b * h(i);
    P(i) = b + 2*h(i);
    R(i) = A(i) / P(i);
end

% Kondisi terhadap waktu -- LF --
Ao = A; An = A;
Qn = Q; Qo = Q;

% Kondisi terhadap waktu -- LW --
Alwo = A; Alwn = A;
Qlwn = Q; Qlwo = Q;


% Simulasi
% Pengulangan waktu
for t = 1:tmax
    
    % Persamaan Kontinuitas
    for i = 2 : imax-1
        An(i) = Ao(i) - (dt/dx)*(Q(i+1)-Q(i-1));
    end

    An(1) = An(2);
    An(imax) = An(imax-1);
    
    % Persamaan Momentum
    for i = 2: imax-1
        R(i) = A(i)/(b+2*(A(i)/b));
        I2 = beta * (((Q(i+1)^2/A(i+1))-(Q(i-1)^2/A(i-1))) / (2*dx));
        I3 = g * A(i) * ((h(i+1)+z(i+1))-(h(i-1)+z(i-1))) / (2*dx);
        I4 = g * abs(Q(i)) * Q(i) * n^2 / (A(i) * R(i)^(4/3));
        Qn(i) = Qo(i) - (I2 + I3 + I4)*2*dt;
    end
    
    Qn(1) = Qn(2);
    Qn(imax) = Qn(imax-1);
   
    % Current -> Initial
    Ao = A;
    Qo = Q;
    % Next -> Current
    A = An;
    Q = Qn;
    
    for i = 1:imax
       h(i) = A(i)/b; 
    end
    
    % Solusi Analitis
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
    
   
    if mod(t, tout) == 0
        p1 = plot(h,'blue');
        hold on
        p2 = plot(hex,'blue');
        hold on
        
        axis manual
        axis ([1 imax 0 2]);
        title('Dam Break Leap-Frog Method');
        xlabel('Panjang(m)');
        ylabel('Elevasi(m)');
        pause(0.001)
        hold off
    end
end