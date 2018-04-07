function object = uma
    object.APR=@APR;
    object.secant=@secant;
    object.pasut=@tidalwave;
    
    object.tvdr=@tvdr;
    object.tvdpsi=@tvdpsi;
    object.tvdcourant=@tvdcourant;
    object.tvdc=@tvdc;
    object.tvdk=@tvdk;
    object.tvdrarr=@tvdrarr;
end



function [A, P, R] = APR(b, m, y)
    A = (b + m*y)*y;
    P = b + 2*y*sqrt(1+m^2);
    R = A/P;
end

function y = secant(yawal, dy, yit, Abm)
A = Abm(1); b = Abm(2); m = Abm(3);
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

function h = tidalwave(t, T, A)
% t = Time
% T = Period
% A = Amplitude
h = sin(2*pi*t/T)*A;
end

function r = tvdr(check, Ar, Qr)
    A = Ar(2); Amin = Ar(1); Aplus = Ar(3);
    Q = Qr(2); Qmin = Qr(1); Qplus = Qr(3);
    
    if check % True for rplus, False for rmin
       r = ((A-Amin)*(Aplus-A)+(Q-Qmin)*(Qplus-Q))/((Aplus-A)*(Aplus-A)+(Qplus-Q)*(Qplus-Q));
    else
       r = ((A-Amin)*(Aplus-A)+(Q-Qmin)*(Qplus-Q))/((A-Amin)*(A-Amin)+(Q-Qmin)*(Q-Qmin));
    end
end

function r = tvdrarr(check, Ar, Qr, i)
    A = Ar(i); Amin = Ar(i-1); Aplus = Ar(i+1);
    Q = Qr(i); Qmin = Qr(i-1); Qplus = Qr(i+1);
    
    if check % True for rplus, False for rmin
       r = ((A-Amin)*(Aplus-A)+(Q-Qmin)*(Qplus-Q))/((Aplus-A)*(Aplus-A)+(Qplus-Q)*(Qplus-Q));
    else
       r = ((A-Amin)*(Aplus-A)+(Q-Qmin)*(Qplus-Q))/((A-Amin)*(A-Amin)+(Q-Qmin)*(Q-Qmin));
    end
end

function psi = tvdpsi(r)
    if r > 0
        psi = min([2*r,1]);
    else
        psi = 0;
    end
end

function courant = tvdcourant(Q, A, B, dt, dx)
    g = 9.81;
    nu = abs(Q)/A + sqrt(g*A/B);
    courant = nu*dt/dx;
end

function c = tvdc(courant)
    if courant <= 0.5
        c = courant*(1-courant);
    else
        c = 0.25;
    end
end

function k = tvdk(r, c)
    k = 0.5*c*(1-tvdpsi(r));
end