%% INFORMATION! READ ME FIRST
% Program Saint-Venant: MacCormack Scheme
% Case: Dambreak
% Script by Taruma (25017046), hi@taruma.info
% Documentation: N/A
% Github this project: https://github.com/taruma/belajar-tsa
% Folder: lika-dambreak/matlab
% WARNING, this script is using Function in script, Use 2016b or later!!

close all;
clc; clear;

uma = uma;

g = 9.81;
n = 0;
beta = 1;

dx = 0.5;
dt = 0.05;

time = 50;
tout = 1;

% Kondisi Dam
b = 1;
h1 = 1.5;
h2 = 1;
lf = 55;
m = 0;

% Grid
imax = lf/dx+1;
tmax = time/dt;

set0 = zeros(1,imax);

Q = set0; z = set0; A = set0; h = set0; R = set0; P = set0;
Qp = set0; zp = set0; Ap = set0; hp = set0; Rp = set0; Pp = set0;
Qc = set0; zc = set0; Ac = set0; hc = set0; Rc = set0; Pc = set0;
Qn = set0; An = set0;

yawal = 5; dy=0.1; yit=10;

%kondisi awal
for i=1:imax
    if i <= imax/2
        h(i) = h1;
    else
        h(i) = h2;
    end
    [A(i), P(i), R(i)] = uma.APR(b, m, h(i));
end

for t = 1:tmax
    %predictor step
    
    %initial
    for i=1:imax
        [A(i), P(i), R(i)] = uma.APR(b, m, h(i));
    end
    
    for i = 2:imax-1
        % Kontinuitas
        Ap(i) = A(i) - dt/dx*(Q(i+1)-Q(i));
        
        % Momentum
        I2p = beta/dx * ((Q(i+1)^2/A(i+1)) - (Q(i)^2/A(i)));
        I3p = g * A(i)/dx * ((h(i+1)+z(i+1))-(h(i)+z(i)));
        I4p = g * Q(i) * abs(Q(i)) * n^2 / (A(i)*R(i)^4/3);
        
        Qp(i) = Q(i) - dt*(I2p+I3p+I4p);
    end
    
    % Syarat batas predictor
    Qp(1) = Qp(2);
    Qp(imax) = Qp(imax-1);
    Ap(1) = Ap(2);
    Ap(imax) = Ap(imax-1);
    
    zp = z;
    
    for i = 1:imax
        hp(i) = uma.secant(yawal, dy, yit, [Ap(i),b,m]);
        [Ap(i),Pp(i),Rp(i)] = uma.APR(b, m, hp(i));
    end
    
    %corrector step
    for i = 2:imax-1
        % Kontinuitas
        Ac(i) = A(i) - dt/dx*(Qp(i)-Qp(i-1));
        
        % Momentum
        I2c = beta/dx * ((Qp(i)^2/Ap(i)) - (Qp(i-1)^2/Ap(i-1)));
        I3c = g * Ap(i) * ((hp(i)+zp(i))-(hp(i-1)+zp(i-1)))/dx;
        I4c = g * Qp(i) * abs(Qp(i)) * n^2 / (Ap(i)*Rp(i)^4/3);
        
        Qc(i) = (Q(i)+Qp(i))/2 - dt*(I2c+I3c+I4c);
    end
    
    % Syarat batas corrector
    Qc(1) = Qc(2);
    Qc(imax) = Qc(imax-1);
    Ac(1) = Ac(2);
    Ac(imax) = Ac(imax-1);
    
    zc = zp;
    
    for i = 1:imax
        hc(i) = uma.secant(yawal, dy, yit, [Ac(i), b, m]);
        [Ac(i), Pc(i), Rc(i)] = uma.APR(b, m, hc(i));
    end
    
    An = Ac;
    Qn = Qc;
    
    for i = 3:imax-2
        courant = uma.tvdcourant(Q(i), A(i), b, dt, dx);
        c = uma.tvdc(courant);
        
        rplus_0 = uma.tvdr(true , [A(i-1),A(i),A(i+1)], [Q(i-1), Q(i), Q(i+1)]);
        rmin_1  = uma.tvdr(false, [A(i-1+1),A(i+1),A(i+1+1)], [Q(i-1+1), Q(i+1), Q(i+1+1)]);
        rplus_1 = uma.tvdr(true , [A(i-2),A(i-1),A(i)], [Q(i-2), Q(i-1), Q(i)]);
        rmin_0  = uma.tvdr(false, [A(i-1),A(i),A(i+1)], [Q(i-1), Q(i), Q(i+1)]);
        
        
        Kp_0 = uma.tvdk(rplus_0, c);
        Km_0 = uma.tvdk(rmin_1, c);
        Kp_1 = uma.tvdk(rplus_1, c);
        Km_1 = uma.tvdk(rmin_0, c);
        
        TVD_A = (Kp_0 + Km_0)*(A(i+1)-A(i))-(Kp_1 + Km_1)*(A(i)-A(i-1));
        TVD_Q = (Kp_0 + Km_0)*(Q(i+1)-Q(i))-(Kp_1 + Km_1)*(Q(i)-Q(i-1));
        
        An(i) = An(i) + TVD_A;
        Qn(i) = An(i) + TVD_Q;
    end
    
    %Syaratbatas
%     An(1) = An(2);
%     An(imax) = An(imax-1);
%     Qn(1) = Qn(2);
%     Qn(imax) = Qn(imax-1);
%     
    A = An;
    Q = Qn;
    
    for i=1:imax
        h(i) = uma.secant(yawal, dy, yit, [A(i), b, m]);
    end
    
        if mod(t, tout) == 0
        p1 = plot(h,'red');
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