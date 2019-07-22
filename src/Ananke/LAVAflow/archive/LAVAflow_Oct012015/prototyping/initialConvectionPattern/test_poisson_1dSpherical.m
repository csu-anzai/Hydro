clear all; close all; clc;

nX = 100;


% =============================================
%              Domain information
% =============================================

xbnds = [1e0,2e0];
xbasis = linspace(xbnds(1),xbnds(2),nX+4);


% =============================================
%              Domain information
% =============================================

potential = zeros(nX+4,1);
RHS = zeros(nX+4,1);
v = zeros(nX+4,1);




% =============================================
%              Initial information
% =============================================
% Initial velocity
% v(:) = 1./(xbasis.^2+0.1);
v(:) = 1./(xbasis);% + 1./(xbasis.^2);
% v(:) = (xbasis);
% v(:) = 1.0;

% Initial divergence
for i=2:nX+3        
    betaX = xbasis(i+1) - xbasis(i-1);
    r    = xbasis(i);
    rSqr = r;
    rP1  = xbasis(i+1);
    rM1  = xbasis(i-1);
    gradX = (1/rSqr)*(rP1*rP1*v(i+1) - rM1*rM1*v(i-1))/betaX;

    RHS(i) = gradX;
end

% potential(:) = rand(nX+4,1);

% RHS(1) = 1;
% RHS(end) = 0.5;

iterMax = 1000000;
errorThresh = 1e-7;
sorCoeff = 1.25;
errorCoords = -1;

minOld = 0;

gradStart = -1.0;
gradEnd = 100.0;



figure
hold on;

% Iteration Loop
for iter=1:iterMax
    
    errorMax = -1;
    
    
%     
%     % Set the lower boundary to open
%     potential(2) = 0;
%     potential(1) = 0;

    % Set the upper boundary to open
%     potential(end) = 0;
%     potential(end-1)= 0;
    
    
    
%     % Set the lower boundary to fixed flow
%     potential(2) = 0;
%     potential(1) = potential(3);

%     
%     % Set the upper boundary to fixed flow
%     potential(end-1) = 0;
%     potential(end) = potential(end-2);
    
%     % Symmetric BC
%     potential(end)   = potential(end-3);
%     potential(end-1) = potential(end-2);
%     potential(1)     = potential(4);
%     potential(2)     = potential(3);
    
%     potential(end) = -potential(end-1);
    potential(1)   = potential(2);
%     potential(end) = potential(end-1);
%     potential(1) = 0;
%     potential(end) = 0;
    potential(end) = -potential(end-1);
    
    
    if(rem(iter,1)==0 && true)
    % X loop
    for i=2:nX+3
        dx = xbasis(i+1)-xbasis(i);
        dxSqr = dx*dx;

        potLapla = 0.0;
        source = 0.0;

        r = xbasis(i);
        rSqr = r*r;
        rM1 = xbasis(i-1);
        rP1 = xbasis(i+1);
        
        % One stencil
        alpha = dx/r;
        potLapla = ((r+dx)/r*potential(i+1)+(r-dx)/r*potential(i-1));
        denom = -2;
        source = dxSqr*RHS(i);

        % New potential value
        potNew = (source - potLapla)/denom;
        % Get the old values
        potOld = potential(i);
        potNew = (1-sorCoeff)*potOld + sorCoeff*potNew;

        err = abs(potNew-potOld)/abs(potNew);

        if(err>errorMax && ~isnan(err))
            errorMax = err;
            errorCoords = [i];
        end


        potential(i) = potNew;

    end
    end 
    
%     % X loop
%     for i=3:1:nX+2
%         dx = xbasis(i+1)-xbasis(i);
%         dxSqr = dx*dx;
% 
%         potLapla = 0.0;
%         source = 0.0;
% 
%         r = xbasis(i);
%         rSqr = r*r;
%         rM1 = xbasis(i-1);
%         rP1 = xbasis(i+1);
% 
% %         % Original
% %         potLapla = potLapla + ((rP1*rP1 - rM1*rM1)/rSqr)*(potential(i+1)-potential(i-1));
% %         potLapla = potLapla + (potential(i+2)+potential(i-2));
% %         % prefactor on the IJ term
% %         denom = -2;
%                 
%         % New
%         potLapla = potLapla + (rP1*rP1*potential(i+2) + rM1*rM1*potential(i-2))/rSqr;
%         denom = -(rP1*rP1 + rM1*rM1)/rSqr;
%                 
%         % source terms
%         source = source + 4.0*dxSqr*RHS(i);
%         
%         % New potential value
%         potNew = (source - potLapla)/denom;
%         % Get the old values
%         potOld = potential(i);
%         potNew = (1-sorCoeff)*potOld + sorCoeff*potNew;
%         err = abs(potNew-potOld)/abs(potNew);
% 
%         if(err>errorMax && ~isnan(err))
%             errorMax = err;
%             errorCoords = [i];
%         end
% 
%         potential(i) = potNew;
% 
%     end
    
    
%     if(rem(iter,5)==0)
%         potential = smooth(potential,3);
%     end
    
    
    gradPot = zeros(nX+4,1);

    for i=2:nX+3
        dx = xbasis(i+1)-xbasis(i-1);
        gradPot(i) = (potential(i+1)-potential(i-1))/dx;
    end
    betaX = xbasis(3) - xbasis(1);
    gradPot(1) =  -(-potential(3) + 4.0*potential(2) - 3.0*potential(1))/betaX;
    betaX = xbasis(end) - xbasis(end-2);
    gradPot(end) =  (potential(end-2) - 4.0*potential(end-1) + 3.0*potential(end))/betaX;
    cla;
    plot(xbasis(3:nX+2),v(3:nX+2),'g--')
    plot(xbasis(3:nX+2),v(3:nX+2)-gradPot(3:nX+2),'b.')
    plot(xbasis(3:nX+2),1./(xbasis(3:nX+2).^2),'ro')
%     plotyy(xbasis,potential,xbasis,gradPot)
ylim([0,2]);
    drawnow;
    

%     minOld = min(potential);
%     potential = potential - minOld;
    
    
    if(rem(iter,1000)==0 )%|| iter == 1)
        fprintf('iteration %d\n',iter);
        fprintf('errorMax: %3.3E (%d)\n',errorMax,errorCoords);
    end
    
    if(errorMax<errorThresh)
        fprintf('Converged! %d %E (%E)\n',iter, errorMax, errorThresh);
        break;
    end
    
end

% Compute the laplacian of the final answer
laplaFinal = zeros(nX+4,1);
for i=2:nX+3
    r = xbasis(i);
    dx = xbasis(i+1)-xbasis(i);
    dxSqr = dx*dx;

%     laplaFinal(i) = (potential(i+1)-2*potential(i)+potential(i-2))/(4*dxSqr);
%     laplaFinal(i) = laplaFinal(i) + 
%     
    laplaFinal(i) = (potential(i+1) - 2*potential(i) + potential(i-1))/dxSqr;
    laplaFinal(i) = laplaFinal(i) + (2/r)*(potential(i+1)-potential(i-1))/(2*dx);
end


%% Compute the gradient of the potential
gradPot = zeros(nX+4,1);

for i=2:nX+3
    dx = xbasis(i+1)-xbasis(i-1);
    gradPot(i) = (potential(i+1)-potential(i-1))/dx;
end
betaX = xbasis(3) - xbasis(1);
gradPot(1) =  -(-potential(3) + 4.0*potential(2) - 3.0*potential(1))/betaX;

betaX = xbasis(end) - xbasis(end-2);
gradPot(end) =  (potential(end-2) - 4.0*potential(end-1) + 3.0*potential(end))/betaX;

% gradPot(1) = (potential(2)-potential(1))/(xbasis(2)-xbasis(1));
% gradPot(end) = (potential(end)-potential(end-1))/(xbasis(end)-xbasis(end-1));



%% Compute new velocity
vNew = v - gradPot;

divV = zeros(nX+4,1);
% Final divergence
for i=2:nX+3        
    betaX = xbasis(i+1) - xbasis(i-1);
    r    = xbasis(i);
    rSqr = r;
    rP1  = xbasis(i+1);
    rM1  = xbasis(i-1);
    gradX = (1/rSqr)*(rP1*rP1*vNew(i+1) - rM1*rM1*vNew(i-1))/betaX;

    divV(i) = gradX;
end

%%
figure
    plotyy(xbasis,potential,xbasis,gradPot)
%% 
figure
subplot(1,2,1)
hold on;
plot(xbasis,v,'g--')
plot(xbasis,vNew)
plot(xbasis,1./(xbasis.^2),'r')
subplot(1,2,2)
plot(xbasis(3:nX+2),divV(3:nX+2))
            


