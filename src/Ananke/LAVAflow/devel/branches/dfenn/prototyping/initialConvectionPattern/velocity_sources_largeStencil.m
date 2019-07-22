clear all; close all; clc;

nX = 50;
nY = 50;

velx = zeros(nX+4,nY+4);
vely = zeros(nX+4,nY+4);
dens = ones(nX+4,nY+4);

xbnds = [0,1];
ybnds = [0,1];
xbasis = linspace(xbnds(1),xbnds(2),nX+4);
ybasis = linspace(ybnds(1),ybnds(2),nY+4);

hVector = [xbasis(2)-xbasis(1);
           ybasis(2)-ybasis(1)];

[X,Y] = meshgrid(xbasis,ybasis);
X = X';
Y = Y';


if(1)
    offset = 0.3;
    velx(:,:) = 0.0;
    vely(:,:) = 0.0;
    nSources = 100;
    for i=1:nSources
        vortexPos = [xbnds(1)+offset, ybnds(1)+offset] + rand(1,2).*[diff(xbnds)-2*offset,diff(ybnds)-2*offset];
        R = sqrt((X-vortexPos(1)).^2 + (Y-vortexPos(2)).^2);
        T = atan2((Y-vortexPos(2)),(X-vortexPos(1)));

        amp = (-1)^(i-1);
        
        velxSource =  amp*R.*sin(T);
        velySource = -amp*R.*cos(T);

        % Weight with hyperbolic tangent to drive field to zero near the boundaries

        offset1 = 0.0;
        offset2 = 0.05;
        RNew1 =  50*(R-offset1);
        RNew2 =  50*(R-offset2);

        weight = (tanh(RNew1)+1)/2 - (tanh(RNew2)+1)/2;

        velx = velx + velxSource.*weight;
        vely = vely + velySource.*weight;

    end
end

if(0)
    dens(:) = 1.0;
end
if(1)
    % Set density field
    densPos = [0.55,0.55];
    dens = 10*exp(-100*((X-densPos(1)).^2 + (Y-densPos(2)).^2))+2;
end
if(0)
    % Set density field
    densPos = [0.5,0.5];
    Lx = 2.0;
    Ly = 2.0;
    amp = 10.0;
    baseline = 5.0+amp;
    
    dens = amp*sin(2*pi*Lx*X).*sin(2*pi*Ly*Y) + baseline;
end

figure
hold on;
contour(X,Y,dens);
quiver(X,Y,velx,vely,3)
axis equal
drawnow



%% Compute the divergence of the velocity field

divV = divergence2(dens.*velx,dens.*vely,xbasis,ybasis,nY+1,nX+1);

fprintf('Initial divergence:\n');
fprintf('\tMin: %E\n',min(min(divV(3:nX+1,3:nY+1))));
fprintf('\tMax: %E\n',max(max(divV(3:nX+1,3:nY+1))));

%% Compute the potential field for divergence


fprintf('Solving for the potential field...\n');
errorThresh = 1e-6;
iterMax = 10000;
divPot = zeros(nX+4,nY+4);
divPot = solvePoissonLargeStencil(divPot, divV,xbasis,ybasis,nX,nY,'dirichlet',errorThresh,iterMax);
% divPot = solvePoisson(divPot, divV,xbasis,ybasis,nX,nY,'dirichlet',errorThresh,iterMax);
% divPot = solvePoisson(divPot, divV,xbasis,ybasis,nX,nY,'outflow',errorThresh,iterMax);


%% Remove the gradient of the divergence potential from the velocity components


fprintf('Removing compressive component...\n');
[gradPotX, gradPotY] = gradient2(divPot, xbasis, ybasis, nX+1,nY+1);

velx = (dens.*velx - gradPotX)./dens;
vely = (dens.*vely - gradPotY)./dens;



%% Compute the divergence of the velocity field again


divV = divergence2(dens.*velx,dens.*vely,xbasis,ybasis,nY+1,nX+1);

fprintf('Final divergence:\n');
fprintf('\tMin: %E\n',min(min(divV(3:nX+1,3:nY+1))));
fprintf('\tMax: %E\n',max(min(divV(3:nX+1,3:nY+1))));


            
%%
figure
hold on;

V = sqrt(velx.^2+vely.^2);
la = V>=0.1*max(V(:));

contour(X,Y,dens);
quiver(X,Y,velx,vely,3)

axis equal


%%



figure
h=pcolor(X,Y,log10(abs(divV)));
set(h,'edgealpha',0);

            

            
            


