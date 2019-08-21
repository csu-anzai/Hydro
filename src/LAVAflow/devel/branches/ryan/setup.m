clear all;
close all;
clc;

ntheta = 200;
nphi = 200;
thetaBasis = linspace(0,2*pi,ntheta+1);
phiBasis = linspace(0,1*pi,nphi+2); 
rbasis = linspace(1,1,1);

[theta,phi,r] = meshgrid(thetaBasis(2:end), phiBasis(2:end-1),rbasis);
x = r.*sin(phi).*cos(theta);
y = r.*sin(theta).*sin(phi);
z = r.*cos(phi);

%dt = DelaunayTri(x(:),y(:),z(:));
%triplot(dt);
T = convhull(x(:),y(:),z(:));
trimesh(T,x(:),y(:),z(:));
axis equal
p = [x(:),y(:),z(:)];

curvature;
meancurv = mean(verteigs);
trisurf(T,p(:,1),p(:,2),p(:,3),meancurv')
colorbar

