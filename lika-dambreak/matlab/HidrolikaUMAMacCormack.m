% Program Saint-venant
% Kasus: Dam Break
% MacCormack Scheme (FD/BD)
% Script by Taruma (25017046), hi@taruma.info
% Dokumentasi: Hidrolika T04 20180222.pdf

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
dt = 0.05;  % time step

time = 15;  % Waktu simulasi
tout = 1;   % Plot tiap waktu tout

% Kondisi Dam
b = 1;      % Lebar saluran
h1 = 1.5;   % Tinggi muka air di kiri (hulu)
h2 = 1;     % Tinggi muka air di kanan (hilir)
lf = 55;    % Panjang saluran

% Grid
imax = lf/dx;   % Jumlah grid untuk jarak (panjang saluran/jarak grid)
tmax = time/dt; % Jumlah time step (waktu simulasi/time step)

% Pendefinisian peubah dan kondisi inisial
% Peubah Corrector
Q = zeros(1, imax); % Debit Q = Q0_0
z = zeros(1, imax); % tinggi datum terhadap dasar permukaan saluran
A = zeros(1, imax); % Luas permukaan
h = zeros(1, imax); % kedalaman air
R = zeros(1, imax); % Jari-jari hidrolis = Luas permukaan / Keliling basah
P = zeros(1, imax); % Keliling basah

% Peubah Predictor (This array for predictor step)
Qp = Q;
zp = z;
Ap = A;
hp = h;

% ------ Solusi Analitik -----
hex = h;            % kedalaman air exact (solusi analitik)

% Kondisi awal
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

% Kondisi terhadap waktu -- LW --
An = A; % An = Ap1_0
Qn = Q; % Qn = Qp1_0

% ------ MacCormack ------ 
% Pengulangan waktu
for t = 1 : tmax
    
    % Predictor Step    
    for i = 2:imax-1
    % Initial
        A(i) = b * h(i);
        P(i) = b + 2*h(i);
        R(i) = A(i) / P(i);
        
        % Kontinuitas
        Ap(i) = A(i) - dt/dx * (Q(i+1)-Q(i));
        
        % Momentum
        I2p = beta * ((Q(i+1)^2/A(i+1)) - (Q(i)^2/A(i)))/dx;
        I3p = g * A(i) * ((h(i+1)+z(i+1))-(h(i)+z(i)))/dx;
        I4p = g * Q(i) * abs(Q(i)) * n^2 / (A(i)*R(i)^4/3);
        Qp(i) = Q(i) - dt * (I2p + I3p + I4p);
    end
    
    % Syarat Batas Predictor
    Qp(1) = Q(2);
    Qp(imax) = Qp(imax-1);
    Ap(1) = Ap(2);
    Ap(imax) = Ap(imax-1);
    
    % Nilai Predictor
    for i = 1:imax
        hp(i) = Ap(i)/b;
        zp(i) = z(i);
    end
    
    % Corrector Step
    for i = 2:imax-1
    % Initial
        A(i) = b * h(i);
        P(i) = b + 2*h(i);
        R(i) = A(i) / P(i);
        
        % Kontinuitas
        Aph = (A(i) + Ap(i))/2;         
        An(i) = Aph - dt/(2*dx) * (Qp(i)-Qp(i-1));
        
        % Momentum
        I2c = beta * ((Qp(i)^2/Ap(i)) - (Qp(i-1)^2/Ap(i-1)))/dx;
        I3c = g * A(i) * ((hp(i)+zp(i))-(hp(i-1)+zp(i-1)))/dx;
        I4c = g * Q(i) * abs(Q(i)) * n^2 / (A(i)*R(i)^4/3);
        
        Qph = (Q(i)+Qp(i))/2;
        Qn(i) = Qph - dt/2 * (I2c + I3c + I4c);
    end    
    
    % Syarat Batas
    Qn(1) = Qn(2);
    Qn(imax) = Qn(imax-1);
    An(1) = An(2);
    An(imax) = An(imax-1);
    
    % Current -> Next Step
    A = An;
    Q = Qn;
    
    % Tinggi H
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