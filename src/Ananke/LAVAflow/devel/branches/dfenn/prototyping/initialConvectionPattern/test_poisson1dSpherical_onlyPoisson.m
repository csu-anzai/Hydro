clear all; close all; clc;

nX = 1000;


% =============================================
%              Domain information
% =============================================

xbnds = [1e0,2e0];
xbasis = linspace(xbnds(1),xbnds(2),nX+4);


% =============================================
%              Domain information
% =============================================

potential = zeros(nX+4,1);
potentialTrue = zeros(nX+4,1);
RHS = zeros(nX+4,1);


% =============================================
%              Initial information
% =============================================
% potential(:) = xbasis;
potentialTrue(:) = xbasis;
potential(:) = (1-2*rand(nX+4,1));
RHS(:) = 2./xbasis;



% =============================================
%              Solve system
% =============================================
iterMax = 100000;
errorThresh = 1e-8;
sorCoeff = 1.75;
errorCoords = -1;

minOld = 0;

gradStart = -1.0;
gradEnd = 100.0;


% Iteration Loop
for iter=1:iterMax
    
    errorMax = -1;
    
    potential(1) = xbasis(1);
    potential(2) = xbasis(1);
    potential(end) = xbasis(end);
    potential(end-1) = xbasis(end-1);
    
    
    
    % X loop
    for i=3:nX+2
        dx = xbasis(i+1)-xbasis(i);
        dxSqr = dx*dx;

        potLapla = 0.0;
        source = 0.0;

        r = xbasis(i);
        rSqr = r*r;
        rM1 = xbasis(i-1);
        rP1 = xbasis(i+1);

        % Original
        potLapla = potLapla + ((rP1*rP1 - rM1*rM1)/rSqr)*(potential(i+1)-potential(i-1));
        potLapla = potLapla + (potential(i+2)+potential(i-2));
        % prefactor on the IJ term
        denom = -2;
        
        
%         % New
%         potLapla = potLapla + (rP1*rP1*potential(i+2) + rM1*rM1*potential(i-2))/rSqr;
%         denom = -(rP1*rP1 + rM1*rM1)/rSqr;
%         
        % source terms
        source = source + 4.0*dxSqr*RHS(i);

%         % One stencil
%         alpha = dx/r;
%         potLapla = ((r+dx)/r*potential(i+1)+(r-dx)/r*potential(i-1));
% %         potLapla = ((1+alpha)*potential(i+1)+(1-alpha)*potential(i-1));
%         denom = -2;
%         source = dxSqr*RHS(i);
        
        % New potential value
        potNew = (source - potLapla)/denom;



        % Get the old values
        potOld = potential(i);
        potNew = (1-sorCoeff)*potOld + sorCoeff*potNew;


        if(~isreal(source))
            fprintf(2,'New potential is complex!\n');
        end

%         err = abs(potNew-potOld-minOld)/abs(potNew-minOld);
        err = abs(potNew-potOld)/abs(potNew);

        if(err>errorMax && ~isnan(err))
            errorMax = err;
            errorCoords = [i];
        end

        potential(i) = potNew;

    end
    
       
    
    if(rem(iter,1000)==0 )%|| iter == 1)
        fprintf('iteration %d\n',iter);
        fprintf('errorMax: %3.3E (%d)\n',errorMax,errorCoords);
    end
    
    if(errorMax<errorThresh)
        fprintf('Converged! %d %E (%E)\n',iter, errorMax, errorThresh);
        break;
    end
    
end

%%

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

figure
hold on;
plot(xbasis(3:nX+2),potential(3:nX+2),'bo')
plot(xbasis(3:nX+2),potentialTrue(3:nX+2),'r-')

