clear all; close all; clc;

nX = 20;
nY = 20;

velr = zeros(nX+4,nY+4);
velp = zeros(nX+4,nY+4);
vorticity = zeros(nX+4,nY+4);



% vorticity(20:40,20:40) = 1.0;


% =============================================
%              Domain information
% =============================================

xbnds = [1e0,1e1];
ybnds = [pi/4,3*pi/4];
xbasis = linspace(xbnds(1),xbnds(2),nX+4);
ybasis = linspace(ybnds(1),ybnds(2),nY+4);

hVector = [xbasis(2)-xbasis(1);
           ybasis(2)-ybasis(1)];

[R,PHI] = meshgrid(xbasis,ybasis);
X = R.*cos(PHI);
Y = R.*sin(PHI);

% =============================================
%              Vorticity sources
% =============================================

if(0)
    vorticity = exp(-100*Rsqr);
end 

% Random Gaussians
if(0)
    offset = 0.3;
    vorticity(:,:) = 0.0;
    nGauss = 100;
    for i=1:nGauss
        pos = [xbnds(1)+offset, ybnds(1)+offset] + rand(1,2).*[diff(xbnds)-2*offset,diff(ybnds)-2*offset];
        Rsqr = (X'-pos(1)).^2 + (Y'-pos(2)).^2;

    %     amplitude = (-1)^i*((posSphere(1)-xbnds(1))/(xbnds(2)-xbnds(1)-offset)*1000 + 0);

        amplitude = (-1)^i;
        sigma = -50;
        vorticity = vorticity + amplitude*exp(sigma*Rsqr);    

    end
end

% Random Spherical Gaussians
if(1)
    charLength = diff(xbnds);
    radOffset = diff(xbnds)/2;%/5;
    phiOffset = .4;
    phiOffset = pi/4;
    nGauss = 1;
    for i=1:nGauss
        posSphere = [xbnds(1)+radOffset, ybnds(1)+phiOffset] + rand(1,2).*[diff(xbnds)-2*radOffset,diff(ybnds)-2*phiOffset];
        pos = [posSphere(1)*cos(posSphere(2)), posSphere(1)*sin(posSphere(2))];
        Rsqr = (X'-pos(1)).^2 + (Y'-pos(2)).^2;

    %     amplitude = (-1)^i*((posSphere(1)-xbnds(1))/(xbnds(2)-xbnds(1)-offset)*1000 + 0);

        amplitude = (-1)^i*1e0;
        sigma = -500;
        
        vorticity = vorticity + amplitude*exp(sigma*Rsqr/charLength^2);
        
    end
end

% Circle of constant vorticity
if(0)
    Rsqr = (X'-0.5).^2 + (Y'-0.5).^2;
    la = Rsqr <= 0.01;
    vorticity(la) = 1.0;
    
%     Rsqr = (X'-0.75).^2 + (Y'-0.75).^2;
%     la = Rsqr <= 0.01;
%     vorticity(la) = -1.0;
end
    


% =============================================
%              Density perturbations
% =============================================
dens = ones(nX+4,nY+4);
gradDensR = zeros(nX+4,nY+4);
gradDensP = zeros(nX+4,nY+4);

if(1)
    dens(:,:) = 1.0;
end

if(0)
    Lr = 1;
    Lp = 1;
    Amp = 1e0;
    dens = Amp*(sin(2*pi*(PHI'-ybnds(1))/(ybnds(2)-ybnds(1))*Lp).*sin(2*pi*(R'-xbnds(1))*Lr/(xbnds(2)-xbnds(1)))+1.01)/2;
end


if(0)
    posSphere = [6e7,pi/2];
    pos = [posSphere(1)*cos(posSphere(2)), posSphere(1)*sin(posSphere(2))];
    Rsqr = (X'-pos(1)).^2 + (Y'-pos(2)).^2;

    amplitude = 1e8;
    sigma = -1000;

    dens = dens + amplitude*exp(sigma*Rsqr/charLength^2);
end



% =============
% Construct grad(dens)
% =============

[gradDensR, gradDensP] = gradient2(dens, xbasis, ybasis, nX+1,nY+1, 'spherical');


% =============================================
%           Solve for velocity field
% =============================================
divergenceSource = zeros(nX+4,nY+4);

iterMax = 10000;

errorMax = 1e99;
errorThresh = 1e-6;

sorCoeff = 0.5;

isForward = true;
% Iteration Loop
for iter=1:iterMax
        
    % Compute the term -dot(u,gradDens)/dens
    divergenceSource = (velr.*gradDensR + velp.*gradDensP)./dens;
    
    % Compute the appropriate source terms required for the velr and velp
    velrSource = zeros(size(velr));
    velpSource = zeros(size(velp));
    for i=2:nX+3
        for j=2:nY+3
            
            betaX = xbasis(i+1) - xbasis(i-1);
            betaY = ybasis(j+1) - ybasis(j-1);
            r = xbasis(i);
            rSqr = r*r;
            phi = ybasis(j);
            
            
%             velrSource(i,j) = velrSource(i,j) + (divergenceSource(i+1,j)-divergenceSource(i-1,j))/betaX;
%             velrSource(i,j) = velrSource(i,j) + -(1/r)*(vorticity(i,j+1)-vorticity(i,j-1))/betaY;
%             velrSource(i,j) = velrSource(i,j) + (2/rSqr)*(velp(i,j+1)-velp(i,j-1))/betaY;
%             velrSource(i,j) = velrSource(i,j) + (cot(phi)/rSqr)*velp(i,j);
%             velrSource(i,j) = velrSource(i,j) + -(cot(phi)/r)*(velp(i+1,j)-velp(i-1,j))/betaX;
%             velrSource(i,j) = velrSource(i,j) + (cot(phi)/r)*(velr(i,j+1)-velr(i,j-1))/betaY;
            
            velrSource(i,j) = velrSource(i,j) + 2/rSqr*(velp(i,j+1)-velp(i,j-1))/betaY;
            velrSource(i,j) = velrSource(i,j) + 2/rSqr*velr(i,j);
            velrSource(i,j) = velrSource(i,j) + cot(phi)/rSqr*(velr(i,j+1)-velr(i,j-1))/betaY;
            velrSource(i,j) = velrSource(i,j) + cot(phi)/rSqr*velp(i,j);
            velrSource(i,j) = velrSource(i,j) - cot(phi)/r*(velp(i+1,j)-velp(i-1,j))/betaX;
            velrSource(i,j) = velrSource(i,j) - (1/r)*(vorticity(i,j+1)-vorticity(i,j-1))/betaY; 
            velrSource(i,j) = velrSource(i,j) - (divergenceSource(i+1,j)-divergenceSource(i-1,j))/betaX; 
%             velrSource(i,j) = velrSource(i,j) +
%             velrSource(i,j) = velrSource(i,j) +
            
            velpSource(i,j) = velpSource(i,j) + (1/r)*(velp(i+1,j)-velp(i-1,j))/betaX;
            velpSource(i,j) = velpSource(i,j) + ((1+csc(phi)^2)/rSqr)*velp(i,j);
            velpSource(i,j) = velpSource(i,j) - (3/rSqr)*(velr(i,j+1)-velr(i,j-1))/betaY;
            
            velpSource(i,j) = velpSource(i,j) - (1/r)*(divergenceSource(i,j+1)-divergenceSource(i,j-1))/betaY;
            velpSource(i,j) = velpSource(i,j) + (vorticity(i+1,j)-vorticity(i-1,j))/betaX;
            
        end
    end
    
    [velr,errorMaxX] = solvePoissonLargeStencil(velr, velrSource,xbasis,ybasis,nX,nY,'dirichlet',errorThresh,1, 'spherical');
    [velp,errorMaxY] = solvePoissonLargeStencil(velp, velpSource,xbasis,ybasis,nX,nY,'dirichlet',errorThresh,1, 'spherical');
%     [velr,errorMaxX] = solvePoissonLargeStencil(velr, velrSource,xbasis,ybasis,nX,nY,'toth',errorThresh,1, 'spherical');
%     [velp,errorMaxY] = solvePoissonLargeStencil(velp, velpSource,xbasis,ybasis,nX,nY,'toth',errorThresh,1, 'spherical');
%     [velr,errorMaxX] = solvePoissonLargeStencil(velr, velrSource,xbasis,ybasis,nX,nY,'outflow',errorThresh,1, 'spherical');
%     [velp,errorMaxY] = solvePoissonLargeStencil(velp, velpSource,xbasis,ybasis,nX,nY,'outflow',errorThresh,1, 'spherical');
%     [velr,errorMaxX] = solvePoissonLargeStencil(velr, velrSource,xbasis,ybasis,nX,nY,'periodic',errorThresh,1, 'spherical');
%     [velp,errorMaxY] = solvePoissonLargeStencil(velp, velpSource,xbasis,ybasis,nX,nY,'periodic',errorThresh,1, 'spherical');
    
    
    if(rem(iter,100)==0 || iter == 1)
        fprintf('iteration %d\n',iter);
        fprintf('errorMaxX: %3.3E\n',errorMaxX);
        fprintf('errorMaxY: %3.3E\n',errorMaxY);
    end
    
    if(errorMaxX<errorThresh && errorMaxY<errorThresh)
        fprintf('Converged! iter = %d\n',iter);
        fprintf('x: %E\ty: %E\tthresh: %E\n',errorMaxX,errorMaxY,errorThresh);
        break;
    end
    
end



% %% Compute the divergence of the velocity field
% % velp = velp./R;
% divV = divergence2(dens.*velr, dens.*velp, xbasis, ybasis, nX+1,nY+1, 'spherical');
% 
% min(divV(:))
% max(divV(:))
% 
% %% Compute the potential field for divergence
% 
% divPot = zeros(nX+4,nY+4);
% % divPot = divV;
% iterMax = 1000;
% errorThresh = 1e-9;
% sorCoeff = 0.75;
% 
% % divPot = solvePoissonLargeStencil(divPot, divV,xbasis,ybasis,nX,nY,'dirichlet',errorThresh,iterMax, 'spherical');
% divPot = solvePoissonLargeStencil(divPot, divV,xbasis,ybasis,nX,nY,'toth',errorThresh,iterMax, 'spherical');
% 
% %% Remove the gradient of the divergence potential from the velocity components
% 
% [gradPotR, gradPotP] = gradient2(divPot, xbasis, ybasis, nX+1,nY+1, 'spherical');
% velr = (dens.*velr - gradPotR)./dens;
% velp = (dens.*velp - gradPotP)./dens;


%% Compute the divergence of the velocity field again

divV = divergence2(dens.*velr, dens.*velp, xbasis, ybasis, nX+1,nY+1, 'spherical');
            
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

% [gradX,gradY] = gradient(pot,hVector(1),hVector(2));
% [X,Y] = meshgrid(xbasis,ybasis);

% [curlZ,~] = curl(X,Y,velx',vely');

V = sqrt(velx.^2+vely.^2);
la = V>=0.01*max(V(:));

contour(X,Y,dens')
quiver(X,Y,velx',vely',1)

axis equal

%%



figure
h=pcolor(X,Y,log10(abs(divV')));
set(h,'edgealpha',0);

min(divV(:))
max(divV(:))
            

            
            


