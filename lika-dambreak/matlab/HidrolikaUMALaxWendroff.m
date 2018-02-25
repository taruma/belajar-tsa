% Program Saint-venant
% Kasus: Dam Break
% Lax Wendroff
% Script by Taruma (25017046), hi@taruma.info
% Dokumentasi: Hidrolika T03 20180220.pdf

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
tout = 0.01;   % Plot tiap waktu tout

% Kondisi Dam
b = 1;      % Lebar saluran
h1 = 1.5;   % Tinggi muka air di kiri (hulu)
h2 = 1;     % Tinggi muka air di kanan (hilir)
lf = 55;    % Panjang saluran

% Grid
imax = lf/dx;   % Jumlah grid untuk jarak (panjang saluran/jarak grid)
tmax = time/dt; % Jumlah time step (waktu simulasi/time step)

% Pendefinisian peubah dan kondisi inisial
% Peubah nol
Q = zeros(1, imax); % Debit Q = Q0_0
z = zeros(1, imax); % tinggi datum terhadap dasar permukaan saluran
A = zeros(1, imax); % Luas permukaan
R = zeros(1, imax); % Jari-jari hidrolis = Luas permukaan / Keliling basah
h = zeros(1, imax); % kedalaman air
P = zeros(1, imax); % Keliling basah

% ------ Solusi Analitik -----
hex = h;            % kedalaman air exact (solusi analitik)

% Kondisi awal
for i = 1:imax
    if i <= imax/2
        h(i) = h1;
    else
        h(i) = h2;
    end
    % A, P, R harus didefinisikan dari awal karena akan digunakan saat mencari nilai A.h_.h dgn . = m/p
    A(i) = b * h(i);    
    P(i) = b + 2*h(i);
    R(i) = A(i) / P(i);
end

% Kondisi terhadap waktu -- LW --
An = A; % An = Ap1_0
Qn = Q; % Qn = Qp1_0

% ------ Lax-Wendroff --------
% Pengulangan waktu
for t = 1 : tmax
    
    for i = 2:imax-1
    % Initial
        A(i) = b * h(i);
        P(i) = b + 2*h(i);
        R(i) = A(i) / P(i);
    
    % Lax Scheme
    % Calculate intermediate mesh points
    % Nilai tengah + - half (1/2)
        A0_ph = (A(i) + A(i+1))/2;
        A0_mh = (A(i) + A(i-1))/2;
        Q0_ph = (Q(i) + Q(i+1))/2;
        Q0_mh = (Q(i) + Q(i-1))/2;
    
    % Kontinuitas
        % Right
        Aph_ph = A0_ph - dt/dx * (Q(i+1)-Q(i))/2;
        % Left
        Aph_mh = A0_mh - dt/dx * (Q(i)-Q(i-1))/2;
    
    % Momentum
        I4 = g * Q(i) * abs(Q(i)) * n^2 / ((A(i) * R(i)^(4/3)));
        % Right
        I2_r = beta * ((Q(i+1)^2/A(i+1)) - (Q(i)^2/A(i)))/dx;
        I3_r = g * A(i) * ((h(i+1)+z(i+1)) - (h(i)+z(i)))/dx; 
        
        Qph_ph = Q0_ph - dt/2 * (I2_r + I3_r + I4);
        
        % Left
        I2_l = beta * ((Q(i)^2/A(i)) - (Q(i-1)^2/A(i-1)))/dx;
        I3_l = g * A(i) * ((h(i)+z(i)) - (h(i-1)+z(i-1)))/dx;
        
        Qph_mh = Q0_mh - dt/2 * (I2_l + I3_l + I4);

    % Nilai i+1/2 dan i-1/2 dan t+1/2
    % Storing value of intermediate mesh points
        hph_ph = Aph_ph/b;
        hph_mh = Aph_mh/b;
        zph_ph = (z(i)+z(i+1))/2;
        zph_mh = (z(i)+z(i-1))/2;
    
    % Wendroff Scheme
    % for i = 2:imax-1
       % Kontinuitas
       An(i) = A(i) - dt/dx *(Qph_ph - Qph_mh);
       
       % Momentum
       I2_w = beta * ((Qph_ph^2/Aph_ph - (Qph_mh^2/Aph_mh)))/dx;
       I3_w = g * A(i) * ((hph_ph+zph_ph) - (hph_mh+zph_mh))/dx;
       
       Qn(i) = Q(i) - dt * (I2_w + I3_w + I4);
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
        title('Dam Break Lax-Wendroff Method');
        xlabel('Grid (x)');
        ylabel('Elevasi(m)');
        pause(0.001)
        hold off
    end
end