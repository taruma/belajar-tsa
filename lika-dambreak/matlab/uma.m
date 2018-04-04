function object = uma
    object.APR=@APR;
    object.secant=@secant;
    object.pasut=@tidalwave;
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
