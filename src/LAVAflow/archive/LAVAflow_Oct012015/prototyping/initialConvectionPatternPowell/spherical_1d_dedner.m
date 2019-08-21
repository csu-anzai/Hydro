clear all; close all; clc;

% ==============================================
%       Grid information
% ==============================================
nxBase = 100;
nGC = 1;
nx = nxBase+2*nGC;
intMin = 1+nGC;
intMax = nGC+nxBase;



% ==============================================
%       Domain information
% ==============================================
xbnds = [1, 10];
lengthX = xbnds(2)-xbnds(1);
h = lengthX/(nx-1);
xbasis = linspace(xbnds(1),xbnds(2),nx);
r = xbasis';
rSqr = r.*r;

% ==============================================
%       Simulation variables
% ==============================================
vx = zeros(nx,1);
vxNew = zeros(nx,1);
divV = zeros(nx,1);
scalar = zeros(nx,1);
scalarNew = zeros(nx,1);
advSpd = zeros(nx,1);

% Evolution variables
cfl = 0.5;
advSpdConst = 0.25*h; % 0.1 cells/sec
advSpd(:) = advSpdConst;
% advSpd(xbasis<0.5) = -advSpd(xbasis<0.5);
dt = cfl*min(h./abs(advSpd));

tStart = 0.0;
tFinal = 100000.123;
time = tStart;
finalPass = false;


% ==============================================
%       Initial field
% ==============================================

% vx(:) = sin(2*pi*4*(xbasis-xbnds(1))/lengthX)+5;
% vx(:) = tanh(10*(xbasis-0.5));
% vx(:) = exp(-10*(xbasis-0.5).^2)+2;
% vx(:) = 1./(rSqr+0.01);
vx(:) = 1./r;
V = vx;
divV(:) = divergence1d(V,xbasis,nx,'spherical');
scalar(:) = divV;
fprintf('%E %E\n',min(log10(abs(divV))),max(log10(abs(divV))));
% scalar = divergence1d(V,xbasis,nx);
% nWavesX = 2.0;
% scalar(:) = sin(2*pi*nWavesX*(x-xbnds(1))/lengthX);
divVold = divV;

figure;
subplot(1,2,1);
subplot(1,2,2);
iter = 0;
% Time loop
while(true)
    iter = iter + 1;
    
    V = vx.*rSqr;
%     divV(:) = divergence1d(V,xbasis,nx)./rSqr;
    divVOld = divergence1d(V,xbasis,nx,'spherical');
    divV = divergence1d(V,xbasis,nx,'spherical');
    
%     if(rem(floor(time),100000)==0)
%         time
%         divV(:) = divergence1d(V,xbasis,nx,'spherical');
%         fprintf('%E %E\n',min(log10(abs(divV))),max(log10(abs(divV))));
%     end
    
    %Draw the plot
    
    subplot(1,2,1);
    [ax,~,~] = plotyy(xbasis,vx,xbasis,scalar);
    set(ax(1),'ylim',[-1,5]);
    %pause(2);
    
    
    subplot(1,2,2);
    plot(xbasis,log10(abs(divVOld)))  
    ylim([-5,0])
    
    drawnow;
    
    
    % Compute the divergence of the velocity field
%     divV(:) = divergence1d(vx,xbasis,nx);
    
    % Evolve field
%     coeff = advSpdConst*dt/h;
    coeff = advSpd.^2.*dt/h;
    
    
    % Compute scalar field development
    scalarNew = scalar;
    scalarNew(intMin:intMax) =  -0.5*coeff(intMin:intMax).*(V(intMin+1:intMax+1)-V(intMin-1:intMax-1)) ...
                                +0.5*coeff(intMin:intMax).*coeff(intMin:intMax).*(V(intMin+1:intMax+1) - 2*V(intMin:intMax) + V(intMin-1:intMax-1));
    scalarNew(intMin:intMax) = scalarNew(intMin:intMax)./rSqr(intMin:intMax);
    scalarNew(intMin:intMax) = scalarNew(intMin:intMax) + scalar(intMin:intMax);
    scalar = scalarNew;
        
        
    % Update velocity field
    if(rem(iter,10)==0)
        vxNew = vx;
        gradScalar = divergence1d(scalar,xbasis,nx,'cartesian');
        vxNew = vxNew - dt*gradScalar;
        vx = vxNew;
    end

    
    % Compute boundary conditions
    % Zero-gradient conditions
%     vx(1) = 1;
%     vx(end) = vx(end-1);
    
    
    % Compute boundary conditions
    % Zero-gradient conditions
%     scalar(1) = 0;
%     scalar(end) = 0;
% scalar(1) = 0;
% scalar(end) = scalar(end-1);
    scalar(1) = scalar(2);
    scalar(end) = scalar(end-1);
    
    
    % Break if this was the last time step
    if(finalPass)
        break;
    end
    
    % Increase simulation time
    timeNew = time + dt;
    if(timeNew>tFinal)
        dt = tFinal-time;
        finalPass = true;
    end
    time = time + dt;
end
    
%%
% Draw the plot
plotyy(xbasis,vx,xbasis,scalar);
ylim([0,4]);

divV(:) = divergence1d(V,xbasis,nx,'spherical');
figure
plot(xbasis,log10(abs(divV)))





