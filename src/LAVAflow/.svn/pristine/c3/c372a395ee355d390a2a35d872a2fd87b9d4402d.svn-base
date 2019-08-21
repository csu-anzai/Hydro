clear all; close all; clc;

nX = 100;
nY = 100;

velx = zeros(nX+2,nY+2);
vely = zeros(nX+2,nY+2);
dens = ones(nX+2,nY+2);

xbnds = [0,1];
ybnds = [0,1];
xbasis = linspace(xbnds(1),xbnds(2),nX+2);
ybasis = linspace(ybnds(1),ybnds(2),nY+2);

hVector = [xbasis(2)-xbasis(1);
           ybasis(2)-ybasis(1)];

[X,Y] = meshgrid(xbasis,ybasis);
X = X';
Y = Y';


if(1)
    offset = 0.5;
    velx(:,:) = 0.0;
    vely(:,:) = 0.0;
    nSources = 1;
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
    % Set density field
    densPos = [0.5,0.5];
    dens = 10*exp(-100*((X-densPos(1)).^2 + (Y-densPos(2)).^2))+2;
end
if(1)
    % Set density field
    densPos = [0.5,0.5];
    Lx = 2.0;
    Ly = 2.0;
    amp = 1.0;
    baseline = 5.0;
    
    dens = amp*sin(2*pi*Lx*X).*sin(2*pi*Ly*Y) + baseline;
end

figure
hold on;
contour(X,Y,dens);
quiver(X,Y,velx,vely,3)
axis equal
drawnow



%% Compute the divergence of the velocity field

divV = divergence2(dens.*velx,dens.*vely,xbasis,ybasis,nY,nX);

fprintf('Initial divergence:\n');
fprintf('\tMin: %E\n',min(divV(:)));
fprintf('\tMax: %E\n',max(divV(:)));

%% Compute the potential field for divergence


fprintf('Solving for the potential field...\n');
errorThresh = 1e-6;
iterMax = 10000;
divPot = zeros(nX+2,nY+2);
divPot = solvePoisson(divPot, divV,xbasis,ybasis,nX,nY,'periodic',errorThresh,iterMax);
% divPot = solvePoisson(divPot, divV,xbasis,ybasis,nX,nY,'dirichlet',errorThresh,iterMax);
% divPot = solvePoisson(divPot, divV,xbasis,ybasis,nX,nY,'outflow',errorThresh,iterMax);


%% Remove the gradient of the divergence potential from the velocity components


fprintf('Removing compressive component...\n');
[gradPotX, gradPotY] = gradient2(divPot, xbasis, ybasis, nX,nY);

velx = (dens.*velx - gradPotX)./dens;
vely = (dens.*vely - gradPotY)./dens;



%% Compute the divergence of the velocity field again


divV = divergence2(dens.*velx,dens.*vely,xbasis,ybasis,nY,nX);

fprintf('Final divergence:\n');
fprintf('\tMin: %E\n',min(divV(:)));
fprintf('\tMax: %E\n',max(divV(:)));


            
%%
figure
hold on;

V = sqrt(velx.^2+vely.^2);
la = V>=0.1*max(V(:));

% contour(X,Y,vorticity')
quiver(X,Y,velx,vely,0)

axis equal


%%



figure
h=pcolor(X,Y,log10(abs(divV)));
set(h,'edgealpha',0);

            

            
            


