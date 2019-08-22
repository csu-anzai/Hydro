clear all; close all; clc;

nX = 100;
nY = 100;

velx = zeros(nX+2,nY+2);
vely = zeros(nX+2,nY+2);

xbnds = [-1,1];
ybnds = [-1,1];
xbasis = linspace(xbnds(1),xbnds(2),nX+2);
ybasis = linspace(ybnds(1),ybnds(2),nY+2);

hVector = [xbasis(2)-xbasis(1);
           ybasis(2)-ybasis(1)];

[X,Y] = meshgrid(xbasis,ybasis);
X = X';
Y = Y';

% Create manufactured solution
RHS = zeros(nX+2,nY+2);
potentialInitial = zeros(nX+2,nY+2);


% ============================
% phi = x^2 + y^2
% => source = 4;
% Dirichlet BC of x^2+y^2
% ============================

RHS(:,:) = 4;
potentialInitial = (X.^2 + Y.^2);
% potentialInitial(2:nX+1,2:nY+1) = 0.0;
potentialAnal = (X.^2 + Y.^2);

% Solve the poisson equation
potential = solvePoisson(potentialInitial, RHS,xbasis,ybasis,nX,nY,'dirichlet', 1e-5,1000);


% Compute error metrics
max(max(sqrt(potential.^2-potentialAnal.^2)))


figure
subplot(2,1,1)
pcolor(X,Y,potential);

subplot(2,1,2)
pcolor(X,Y,potentialAnal);