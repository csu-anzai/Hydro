clear all; close all; clc;

nX = 64;
nY = 64;

velr = zeros(nX+4,nY+4);
velp = zeros(nX+4,nY+4);
velx = zeros(nX+4,nY+4);
vely = zeros(nX+4,nY+4);
    
dens = ones(nX+4,nY+4);


% =============================================
%              Domain information
% =============================================

xbnds = [1e0,2e0];
ybnds = [-pi/3,pi/3] + pi/2;
xbasis = linspace(xbnds(1),xbnds(2),nX+4);
ybasis = linspace(ybnds(1),ybnds(2),nY+4);

hVector = [xbasis(2)-xbasis(1);
           ybasis(2)-ybasis(1)];

[R,PHI] = meshgrid(xbasis,ybasis);
R = R';
PHI = PHI';
X = R.*cos(PHI);
Y = R.*sin(PHI);

% xbnds = [1e9,1.0001e9];
% ybnds = [-pi/5000,pi/5000] + pi/2;
% xbasis = linspace(xbnds(1),xbnds(2),nX+4);
% ybasis = linspace(ybnds(1),ybnds(2),nY+4);
% 
% hVector = [xbasis(2)-xbasis(1);
%            ybasis(2)-ybasis(1)];
% 
% [R,PHI] = meshgrid(xbasis,ybasis);
% R = R';
% PHI = PHI';
% X = R.*cos(PHI);
% Y = R.*sin(PHI);


% =============================================
%              Initial velocity distributions
% =============================================

if(0)

    velr(20:30,20:30) = R(20:30,20:30);
    velp(:,:) = 0.0;

end
if(1)
    rOffset = .4;
    pOffset = pi/8;
    velx(:,:) = 0.0;
    vely(:,:) = 0.0;
    nSources = 100;
    for i=1:nSources
        vortexPosSphere = [xbnds(1)+rOffset, ybnds(1)+pOffset] + rand(1,2).*[diff(xbnds)-2*rOffset,diff(ybnds)-2*pOffset];
        vortexPos = [vortexPosSphere(1)*cos(vortexPosSphere(2)), vortexPosSphere(1)*sin(vortexPosSphere(2))];
        Rtmp = sqrt((X-vortexPos(1)).^2 + (Y-vortexPos(2)).^2);
        T = atan2((Y-vortexPos(2)),(X-vortexPos(1)));

        amp = (-1)^(i-1);
        
        velxSource =  amp*Rtmp.*sin(T);
        velySource = -amp*Rtmp.*cos(T);

        % Weight with hyperbolic tangent to drive field to zero near the boundaries

        offset1 = 0.0;
        offset2 = .2;
        RNew1 = 10*(Rtmp-offset1);
        RNew2 = 10*(Rtmp-offset2);

        weight = (tanh(RNew1)+1)/2 - (tanh(RNew2)+1)/2;

        velx = velx + velxSource.*weight;
        vely = vely + velySource.*weight;
        

    end
    
    
    % X loop
    for i=1:nX+4

        % Y loop
        for j=1:nY+4

            velr(i,j) =  velx(i,j)*cos(ybasis(j)) + vely(i,j)*sin(ybasis(j));
            velp(i,j) = -velx(i,j)*sin(ybasis(j)) + vely(i,j)*cos(ybasis(j));
%             velp(i,j) = velp(i,j)/R(i,j);
        end
    end
end

if(0)
    velr = 1./(R.*R) + 1./R;
    velp(:,:) = 0;
    
    % X loop
    for i=1:nX+4

        % Y loop
        for j=1:nY+4

            velx(i,j) = velr(i,j)*cos(ybasis(j)) - velp(i,j)*sin(ybasis(j));
            vely(i,j) = velr(i,j)*sin(ybasis(j)) + velp(i,j)*cos(ybasis(j));

        end
    end

end


% =============================================
%              Initial density distributions
% =============================================

if(0)
    dens(:,:) = 1.0;
end

if(0)
    Lr = 1.5;
    Lp = 1.5;
    Amp = 1e0;
    dens = Amp*(sin(2*pi*(PHI-ybnds(1))/(ybnds(2)-ybnds(1))*Lp).*sin(2*pi*(R-xbnds(1))*Lr/(xbnds(2)-xbnds(1)))+1.01)/2+5;
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


if(1)
    % Set density field
    Lx = 0.5;
    amp = 10.0;
    baseline = 5.0+amp;
    
    dens = amp*sin(2*pi*Lx*(R-xbnds(1))/(xbnds(2)-xbnds(1))) + baseline;
end


figure
hold on;
% contour(X,Y,dens);
h = pcolor(X,Y,dens); set(h,'edgealpha',0);
quiver(X,Y,velx,vely,1)
axis equal
drawnow



%% Compute the divergence of the velocity field

divV = divergence2(dens.*velr,dens.*velp,xbasis,ybasis,nY+1,nX+1,'spherical');

fprintf('Initial divergence:\n');
fprintf('\tMin: %E\n',min(min(abs(divV(3:nX+2,3:nY+2)))));
fprintf('\tMax: %E\n',max(max(abs(divV(3:nX+2,3:nY+2)))));

%% Compute the potential field for divergence


fprintf('Solving for the potential field...\n');
errorThresh = 1e-3;
iterMax = 10000;
divPot = zeros(nX+4,nY+4);



divPot = solvePoissonLargeStencil(divPot, divV,xbasis,ybasis,nX,nY,'toth',errorThresh,iterMax,'spherical');
% divPot = solvePoissonLargeStencil(divPot, divV,xbasis,ybasis,nX,nY,'periodic',errorThresh,iterMax,'spherical');
% divPot = solvePoissonLargeStencil(divPot, divV,xbasis,ybasis,nX,nY,'dirichlet',errorThresh,iterMax,'spherical');
% divPot = solvePoissonLargeStencil(divPot, divV,xbasis,ybasis,nX,nY,'outflow',errorThresh,iterMax,'spherical');


% if(1)
%     % Set up dot(n,v) conditions for solenoidal boundary conditions
%     % Assume all normals are pointing out of the domain
% 
%     % Bottom
%     divPot(1,1:end) = -velr(1,1:end);
%     divPot(2,1:end) = -velr(2,1:end);
% 
%     % Top
%     divPot(end,1:end)   = velr(end,1:end);
%     divPot(end-1,1:end) = velr(end-1,1:end);
% 
%     % Right
%     divPot(1:end,1) = -velp(1:end,1)./R(1:end,1);
%     divPot(1:end,2) = -velp(1:end,2)./R(1:end,2);
% 
%     % Right
%     divPot(1:end,end)   = velp(1:end,end)./R(1:end,end);
%     divPot(1:end,end-1) = velp(1:end,end-1)./R(1:end,end-1);
% 
%     divPotInitial = divPot;
%     
% figure
% h=pcolor(X,Y,divPot);
% set(h,'edgealpha',0);
% drawnow
% 
%     divPot = solvePoissonLargeStencil(divPot, divV,xbasis,ybasis,nX,nY,'solenoidal',errorThresh,iterMax,'spherical');
%     
% end

%% Remove the gradient of the divergence potential from the velocity components


fprintf('Removing compressive component...\n');
[gradPotX, gradPotY] = gradient2(divPot, xbasis, ybasis, nX+1,nY+1,'spherical');

interiorX = 3:nX+2;
interiorY = 3:nY+2;
velr = (dens.*velr - gradPotX)./dens;
velp = (dens.*velp - gradPotY)./dens;
% velr(interiorX,interiorY) = (dens(interiorX,interiorY).*velr(interiorX,interiorY) - gradPotX(interiorX,interiorY))./dens(interiorX,interiorY);
% velp(interiorX,interiorY) = (dens(interiorX,interiorY).*velp(interiorX,interiorY) - gradPotY(interiorX,interiorY))./dens(interiorX,interiorY);



%% Compute the divergence of the velocity field again


divV = divergence2(dens.*velr,dens.*velp,xbasis,ybasis,nY+1,nX+1,'spherical');

fprintf('Final divergence:\n');
fprintf('\tMin: %E\n',min(min(divV(3:nX+2,3:nY+2))));
fprintf('\tMax: %E\n',max(min(divV(3:nX+2,3:nY+2))));

%% Convert vectors into cartesian
velx = zeros(nX+4,nY+4);
vely = zeros(nX+4,nY+4);

% X loop
for i=1:nX+4

    % Y loop
    for j=1:nY+4
        
        velx(i,j) = velr(i,j)*cos(ybasis(j)) - velp(i,j)*sin(ybasis(j));
        vely(i,j) = velr(i,j)*sin(ybasis(j)) + velp(i,j)*cos(ybasis(j));
        
    end
end

            
%%
figure
hold on;

V = sqrt(velx.^2+vely.^2);
% la = V>=0.1*max(V(:));

contour(X,Y,dens);
quiver(X,Y,velx,vely,1)

axis equal


%%



figure
h=pcolor(X(3:nX+2,3:nY+2),Y(3:nX+2,3:nY+2),log10(abs(divV(3:nX+2,3:nY+2))));
% h=pcolor(X,Y,divPot);
set(h,'edgealpha',0);

            

            
            


