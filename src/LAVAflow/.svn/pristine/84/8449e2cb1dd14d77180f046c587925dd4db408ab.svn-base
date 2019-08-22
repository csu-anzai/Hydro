clear all; close all; clc;

nX = 50;
nY = 50;

velr = zeros(nX+4,nY+4);
velp = zeros(nX+4,nY+4);
velx = zeros(nX+4,nY+4);
vely = zeros(nX+4,nY+4);
    
dens = ones(nX+4,nY+4);


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



if(0)

    velr(20:30,20:30) = R(20:30,20:30);
    velp(:,:) = 0.0;

end
if(0)
    rOffset = 4.5;
    pOffset = pi/4;
    velx(:,:) = 0.0;
    vely(:,:) = 0.0;
    nSources = 1;
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
        offset2 = 2;
        RNew1 = 5*(Rtmp-offset1);
        RNew2 = 5*(Rtmp-offset2);

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
    velr = R;
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

if(1)
    velr = R;
    velp = R.*PHI;
    
    % X loop
    for i=1:nX+4

        % Y loop
        for j=1:nY+4

            velx(i,j) = velr(i,j)*cos(ybasis(j)) - velp(i,j)*sin(ybasis(j));
            vely(i,j) = velr(i,j)*sin(ybasis(j)) + velp(i,j)*cos(ybasis(j));

        end
    end

end

if(1)
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

figure
hold on;
contour(X,Y,dens);
quiver(X,Y,velx,vely,1)
axis equal
drawnow



%% Compute the divergence of the velocity field

divV = divergence2(dens.*velr,dens.*velp,xbasis,ybasis,nY+1,nX+1,'spherical');

fprintf('Initial divergence:\n');
fprintf('\tMin: %E\n',min(min(divV(3:nX+1,3:nY+1))));
fprintf('\tMax: %E\n',max(max(divV(3:nX+1,3:nY+1))));

%% Compare the analyitical divergence to the computed divergence
divAnal = zeros(size(divV));
divAnal = 4+cot(PHI).*PHI;

divDiff = abs(divAnal-divV);
max(max(divDiff(3:nX+2,3:nY+2)))
min(min(divDiff(3:nX+2,3:nY+2)))

