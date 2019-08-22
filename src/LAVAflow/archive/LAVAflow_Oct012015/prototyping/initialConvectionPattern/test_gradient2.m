clear all; close all; clc;

nX = 10;
nY = 10;

scalar = zeros(nX+2,nY+2);

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

% ============================
% phi = x^2 + y^2
% => grad(phi) = [2x, 2y]
% ============================

scalar(:,:) = X.^2 + Y.^2;
gradAnalX = 2*X;
gradAnalY = 2*Y;


% ============================
% phi = x + y
% => grad(phi) = [1, 1]
% ============================

% scalar(:,:) = X + Y;
% gradAnalX = ones(size(scalar));
% gradAnalY = ones(size(scalar));


% Compute gradient
[gradX,gradY] = gradient2(scalar, xbasis, ybasis, nX,nY);

% Compute error metrics

errX = max(abs(gradX(:)-gradAnalX(:)))
errY = max(abs(gradY(:)-gradAnalY(:)))


figure

subplot(2,2,1)
stem3(X,Y,gradX);

subplot(2,2,2)
stem3(X,Y,gradY);


subplot(2,2,3)
pcolor(X,Y,gradAnalX);

subplot(2,2,4)
pcolor(X,Y,gradAnalY);

% Compute the divergence of the vector field
% For phi=x+y the divergence should be 0
% For phi=x^2+y^2 the divergence should be 4
div = divergence2(gradX,gradY,xbasis,ybasis,nY,nX);
min(div(:))
max(div(:))

