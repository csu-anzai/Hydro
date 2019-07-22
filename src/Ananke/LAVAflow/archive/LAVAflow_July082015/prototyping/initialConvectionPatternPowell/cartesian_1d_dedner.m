clear all; close all; clc;

% ==============================================
%       Grid information
% ==============================================
nxBase = 1000;
nGC = 1;
nx = nxBase+2*nGC;
intMin = 1+nGC;
intMax = nGC+nxBase;



% ==============================================
%       Domain information
% ==============================================
xbnds = [0, 1];
lengthX = xbnds(2)-xbnds(1);
h = lengthX/(nx-1);
xbasis = linspace(xbnds(1),xbnds(2),nx);

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
cfl = 0.9;
advSpdConst = 0.01*h; % 0.1 cells/sec
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
vx(:) = exp(-10*(xbasis-0.5).^2)+2;
divV(:) = divergence1d(vx,xbasis,nx)';
% scalar = divergence1d(vx,xbasis,nx);
% nWavesX = 2.0;
% scalar(:) = sin(2*pi*nWavesX*(x-xbnds(1))/lengthX);


figure;

% Time loop
while(true)
    
    
    % Draw the plot
    plotyy(xbasis,vx,xbasis,scalar);
    ylim([0,4]);
%     pause(0.001);
drawnow;
    
    
    % Compute the divergence of the velocity field
%     divV(:) = divergence1d(vx,xbasis,nx);
    
    % Evolve field
    coeff = advSpdConst*dt/h;
    coeff = advSpd.^2.*dt/h;
    
    % Compute scalar field development
    scalarNew = scalar;
    scalarNew(intMin:intMax) =  scalar(intMin:intMax) - ...
                                0.5*coeff(intMin:intMax).*(vx(intMin+1:intMax+1)-vx(intMin-1:intMax-1)) + ...
                                0.5*coeff(intMin:intMax).*coeff(intMin:intMax).*(vx(intMin+1:intMax+1) - 2*vx(intMin:intMax) + vx(intMin-1:intMax-1));
    
    scalar = scalarNew;
    
    
        
    % Update velocity field
    vxNew = vx;
    gradScalar = divergence1d(scalar,xbasis,nx);
%     size(vxNew)
%     size(gradScalar)
    vxNew = vxNew - dt*gradScalar;
    vx = vxNew;

    % Compute boundary conditions
    % Zero-gradient conditions
%     scalar(1) = 0;
%     scalar(end) = 0;
    scalar(1) = scalar(2);
    scalar(end) = scalar(end-1);
    
    % Compute boundary conditions
    % Zero-gradient conditions
%     vx(1) = vx(2);
%     vx(end) = vx(end-1);
    
    
    
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
    
    
% Draw the plot
plotyy(xbasis,vx,xbasis,scalar);
ylim([0,4]);





