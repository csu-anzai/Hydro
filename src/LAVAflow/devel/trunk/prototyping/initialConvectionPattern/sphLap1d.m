clear all; close all; clc;

gradStart = 1;
gradFinal = 10.0;

yStart = 0;
yEnd = 5;

sphLap1D = @(x,f) [f(1); 3-2*f(1)/x];
res = @(y0,y1) [y0(2)-gradStart; y1(2)-gradFinal];

sol = bvpinit(linspace(1,10,100),[10,10]);
sol = bvp4c(sphLap1D,res,sol);

x = linspace(1,10,100);
Sy = deval(sol,x);

% figure
% plotyy(sol.x,sol.y(1,:),sol.x,sol.yp(2,:))

% figure
% plot(x,y)

figure
plotyy(x,Sy(1,:),x,Sy(2,:))