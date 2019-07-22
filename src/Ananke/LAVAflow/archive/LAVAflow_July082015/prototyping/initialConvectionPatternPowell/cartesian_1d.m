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
xbnds = [0, 1];
lengthX = xbnds(2)-xbnds(1);
h = lengthX/(nx-1);
xbasis = linspace(xbnds(1),xbnds(2),nx);

% ==============================================
%       Simulation variables
% ==============================================
vx = zeros(nx,1);
vxNew = zeros(nx,1);
scalar = zeros(nx,1);
scalarNew = zeros(nx,1);
advSpd = zeros(nx,1);

% Evolution variables
cfl = 0.5;
advSpdConst = 0.1*h; % 0.1 cells/sec
advSpd(:) = advSpdConst;
% advSpd(xbasis<0.5) = -advSpd(xbasis<0.5);
dt = cfl*min(h./abs(advSpd));

tStart = 0.0;
tFinal = 1000.123;
time = tStart;
finalPass = false;


% ==============================================
%       Initial field
% ==============================================

% vx = sin(2*pi*4*(xbasis-xbnds(1))/lengthX)+5;
vx = tanh(10*(xbasis-0.5));
scalar = divergence1d(vx,xbasis,nx);
% nWavesX = 2.0;
% scalar(:) = sin(2*pi*nWavesX*(x-xbnds(1))/lengthX);


figure;

% Time loop
while(true)
    
    
    % Evolve field
    coeff = advSpdConst*dt/h;
    coeff = advSpd*dt/h;
    
%     scalarNew = scalar;
%     scalarNew(intMin:intMax) =  scalar(intMin:intMax) - ...
%                                 0.5*coeff(intMin:intMax)'.*(scalar(intMin+1:intMax+1)-scalar(intMin-1:intMax-1)) + ...
%                                 0.5*coeff(intMin:intMax)'.*coeff(intMin:intMax)'.*(scalar(intMin+1:intMax+1) - 2*scalar(intMin:intMax) + scalar(intMin-1:intMax-1));
%     
%     scalar = scalarNew;
%     
%     
%     % Compute boundary conditions
%     % Zero-gradient conditions
%     scalar(1) = scalar(2);
%     scalar(end) = scalar(end-1);
%     
%     % Zero dirichlet conditions
%     scalar(1) = 0;
%     scalar(end) = 0;
    
    
    vxNew = vx;
    vxNew(intMin:intMax) =  vx(intMin:intMax) - ...
                            0.5*coeff(intMin:intMax)'.*(vx(intMin+1:intMax+1)-vx(intMin-1:intMax-1)) + ...
                            0.5*coeff(intMin:intMax)'.*coeff(intMin:intMax)'.*(vx(intMin+1:intMax+1) - 2*vx(intMin:intMax) + vx(intMin-1:intMax-1)) - ...
                            dt.*scalar(intMin:intMax)*advSpd(intMin:intMax);
    
    vx = vxNew;
    
    scalar = divergence1d(vx,xbasis,nx);
    
    % Compute boundary conditions
    % Zero-gradient conditions
    vx(1) = vx(2);
    vx(end) = vx(end-1);
    scalar(1) = scalar(2);
    scalar(end) = scalar(end-1);
    
    % Zero dirichlet conditions
%     vx(1) = 0;
%     vx(end) = 0;
    
    % Draw the plot
    plotyy(xbasis,vx,xbasis,scalar);
    pause(0.05);
    
    
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
    
    




