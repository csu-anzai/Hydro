clc, clear all, clf

% function to test
f = @(x,y,z) 12.7*sin(x) + pi*y.^3 - z.^3;

% integral of function
I = @(xmin,xmax,ymin,ymax,zmin,zmax) triplequad(f,xmin,xmax,ymin,ymax,zmin,zmax,1e-15);

% grid params
xmin = 0; xmax = 1;
ymin = 0; ymax = 1;
zmin = 0; zmax = 1;

% create vectors of spatial data
xc = linspace(xmin,xmax,4);
yc = linspace(ymin,ymax,4);
zc = linspace(zmin,zmax,4);
dx = abs(xc(2) - xc(1));

% initialze cube of data (the interpolation stencil)
uc = zeros(3,3,3);

% fill data with cell averages
for i = 1:length(xc)-1
    for j = 1:length(yc)-1
        for k = 1:length(zc)-1
            uc(i,j,k) = I(xc(i),xc(i+1),yc(i),yc(i+1),zc(i),zc(i+1)) / dx^3;
        end
    end
end

% compute interpolation and detail coefficients
d1 = abs(uc(2,2,2) - 1/8*(uc(3,2,2)-uc(1,2,2)) - 1/8*(uc(2,3,2)-uc(2,1,2)) ...
     - 1/8*(uc(2,2,3)-uc(2,2,1)) + 1/64*(uc(3,3,2) - uc(1,3,2) + uc(3,2,3) ...
     - uc(1,2,3) + uc(2,3,3) - uc(2,1,3) + uc(2,1,1) - uc(2,3,1) + uc(1,2,1) ...
     - uc(3,2,1) + uc(1,1,2) - uc(3,1,2)) + 1/512*(uc(3,3,1) - uc(3,3,3) ...
     + uc(3,1,3) - uc(3,1,1) + uc(1,3,3) - uc(1,3,1) + uc(1,1,1) - uc(1,1,3)) ...
     - I(xc(2),0.5*(xc(2)+xc(3)),yc(2),0.5*(yc(2)+yc(3)),zc(2),0.5*(zc(2)+zc(3)))/(dx/2)^3)

d2 = abs(uc(2,2,2) - 1/8*(uc(3,2,2)-uc(1,2,2)) - 1/8*(uc(2,3,2)-uc(2,1,2)) ...
     - 1/8*(uc(2,2,3)-uc(2,2,1)) + 1/64*(uc(3,3,2) - uc(1,3,2) + uc(3,2,3) ...
     - uc(1,2,3) + uc(2,3,3) - uc(2,1,3) + uc(2,1,1) - uc(2,3,1) + uc(1,2,1) ...
     - uc(3,2,1) + uc(1,1,2) - uc(3,1,2)) + 1/512*(uc(3,3,1) - uc(3,3,3) ...
     + uc(3,1,3) - uc(3,1,1) + uc(1,3,3) - uc(1,3,1) + uc(1,1,1) - uc(1,1,3)) ...
     - I(0.5*(xc(2)+xc(3)),xc(3),yc(2),0.5*(yc(2)+yc(3)),zc(2),0.5*(zc(2)+zc(3)))/(dx/2)^3)
