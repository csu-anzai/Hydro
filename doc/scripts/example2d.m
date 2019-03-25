clc, clear all, clf

% number of cells
n = 32;

% create mesh
h = 1 / n;
xh = -h:h:1+h
yh = -h:h:1+h
xh = xh(1:end-1) + h/2;
yh = yh(1:end-1) + h/2;
h = 1 / (n/2);
x2h = -h:h:1+h
y2h = -h:h:1+h
x2h = x2h(1:end-1) + h/2;
y2h = y2h(1:end-1) + h/2;

% function and integral
f = @(x,y) exp(-4*(x^2 + y^2)) * sin(pi*x) * sin(pi*y);

% wavelet decomposition of function
for i = 1 : n/2
    for j = 1 : n/2
        uij = f( x2h(i), y2h(j) );
        u1 = f( x2h(i)+h, y2h(j) );
        u2 = f( x2h(i)-h, y2h(j) );
        u3 = f( x2h(i), y2h(j)+h );
        u4 = f( x2h(i), y2h(j)-h );
        u5 = f( x2h(i)+h, y2h(j)+h );
        u6 = f( x2h(i)+h, y2h(j)-h );
        u7 = f( x2h(i)-h, y2h(j)+h );
        u8 = f( x2h(i)-h, y2h(j)-h );
        detail(i,j) = f(x2h(i)+h/2,y2h(j)+h/2) - (uij - (u1-u2)/8 - (u3-u4)/8 ...
                        - ( (u5-u6) - (u7-u8) )/64 );
    end
end

[X,Y] = meshgrid( linspace(0,1,16), linspace(0,1,16) );
figure(1)
contourf( X, Y, abs(detail) );
