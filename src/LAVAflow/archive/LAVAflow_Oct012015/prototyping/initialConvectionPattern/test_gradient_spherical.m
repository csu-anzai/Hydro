clear all; close all; clc;

nX = 50;
nY = 50;

    
scalar = zeros(nX+4,nY+4);


% =============================================
%              Domain information
% =============================================

xbnds = [1e0,10];
ybnds = [pi/4,3*pi/4];
xbasis = linspace(xbnds(1),xbnds(2),nX+4);
ybasis = linspace(ybnds(1),ybnds(2),nY+4);

hVector = [xbasis(2)-xbasis(1);
           ybasis(2)-ybasis(1)];

[R,PHI] = meshgrid(xbasis,ybasis);
R = R';
PHI = PHI';
X = R.*cos(PHI);
Y = R.*sin(PHI);


scalar = R + PHI;



%% Compute the gradient of the scalar field

[gradR,gradP] = gradient2(scalar, xbasis, ybasis, nX+1,nY+1, 'spherical');


%% Compare the analyitical divergence to the computed divergence
gradRAnal = ones(size(scalar));
gradPAnal = 1./R;


divDiff = abs(gradRAnal-gradR);
max(max(divDiff(3:nX+2,3:nY+2)))
min(min(divDiff(3:nX+2,3:nY+2)))

divDiff = abs(gradPAnal-gradP);
max(max(divDiff(3:nX+2,3:nY+2)))
min(min(divDiff(3:nX+2,3:nY+2)))

