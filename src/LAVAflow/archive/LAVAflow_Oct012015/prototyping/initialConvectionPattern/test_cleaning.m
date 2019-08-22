clear all; close all; clc;



nX = 101;
nY = 101;

velx = zeros(nX+2,nY+2);
vely = zeros(nX+2,nY+2);

xbnds = [-1,1];
ybnds = [-1,1];
xbasis = linspace(xbnds(1),xbnds(2),nX+2);
ybasis = linspace(ybnds(1),ybnds(2),nY+2);

hVector = [xbasis(2)-xbasis(1);
           ybasis(2)-ybasis(1)];

[X,Y] = meshgrid(xbasis,ybasis);
Rsqr = (X'-0.0).^2 + (Y'-0.0).^2;

R = sqrt(X.^2 + Y.^2)';
T = atan2(Y,X);

velx = R'.*sin(T') + cos(T');
vely = -R'.*cos(T') + sin(T');

% velx = (1./R').*sin(T');% + cos(T');
% vely = -(1./R').*cos(T');% + sin(T');


% Weight with hyperbolic tangent to drive field to zero near the boundaries

offset1 = 0.0;
offset2 = 10;
RNew1 =  25*(R-offset1);
RNew2 =  25*(R-offset2);

weight = (tanh(RNew1)+1)/2 - (tanh(RNew2)+1)/2;

velx = velx.*weight;
vely = vely.*weight;

% velx = cos(T');
% vely = sin(T');

% velx = X'+1;
% vely(:,:) = 0;


% velx = sin(2*pi*X');
% vely(:,:) = 0;

figure
hold on;

quiver(X,Y,velx',vely',1)

% velx(1,:) = 0;
% velx(end,:) = 0;
% velx(:,1) = 0;
% velx(:,end) = 0;
% vely(1,:) = 0;
% vely(end,:) = 0;
% vely(:,1) = 0;
% vely(:,end) = 0;


%% Compute the divergence of the velocity field

divV = zeros(nX+2,nY+2);

% X loop
for i=2:nX+1

    % Y loop
    for j=2:nY+1
        
        % http://ocw.usu.edu/civil_and_environmental_engineering/numerical_methods_in_civil_engineering/LaplaceDirichletTemperature.pdf
        
        % Compute the gradients
        betaX = xbasis(i+1) - xbasis(i-1);
        betaY = ybasis(j+1) - ybasis(j-1);
        
        % Compute the source terms
        gradX =  (velx(i+1,j) - velx(i-1,j))/betaX;
        gradY =  (vely(i,j+1) - vely(i,j-1))/betaY;
        
        divV(i,j) = gradX + gradY;

    end
end

% Compute the divergence on the boundary via forward/backward differences

% Left boundary
for i=1:1
    for j=2:nY+1
        
        % Compute the gradients
        betaX = xbasis(i+2) - xbasis(i);
        betaY = ybasis(j+1) - ybasis(j-1);
        
        % Compute the source terms
        gradX =  (-velx(i+2,j) + 4.0*velx(i+1,j) - 3.0*velx(i,j))/betaX;
        gradY =  (vely(i,j+1) - vely(i,j-1))/betaY;
        divV(i,j) = gradX + gradY;
    end
end


% Right boundary
for i=nX+2:nX+2
    for j=2:nY+1
        
        % Compute the gradients
        betaX = xbasis(i) - xbasis(i-2);
        betaY = ybasis(j+1) - ybasis(j-1);
        
        % Compute the source terms
        gradX =  (velx(i-2,j) - 4.0*velx(i-1,j) + 3.0*velx(i,j))/betaX;
        gradY =  (vely(i,j+1) - vely(i,j-1))/betaY;
        divV(i,j) = gradX + gradY;
    end
end


% Bottom boundary
for i=2:nX+1
    for j=1:1
        
        % Compute the gradients
        betaX = xbasis(i+1) - xbasis(i-1);
        betaY = ybasis(j+2) - ybasis(j);
        
        % Compute the source terms
        gradX =  (velx(i+1,j) - velx(i-1,j))/betaX;
        gradY =  (-vely(i,j+2) + 4.0*vely(i,j+1) - 3.0*vely(i,j))/betaY;
        divV(i,j) = gradX + gradY;
    end
end


% Top boundary
for i=2:nX+1
    for j=nY+2:nY+2
        
        % Compute the gradients
        betaX = xbasis(i+1) - xbasis(i-1);
        betaY = ybasis(j) - ybasis(j-2);
        
        % Compute the source terms
        gradX =  (velx(i+1,j) - velx(i-1,j))/betaX;
        gradY =  (vely(i,j-2) - 4.0*vely(i,j-1) + 3.0*vely(i,j))/betaY;
        divV(i,j) = gradX + gradY;
    end
end

min(divV(:))
max(divV(:))

V = sqrt(velx.^2+vely.^2);
max(V(:))

divVanal = 2*pi*cos(2*pi*X');

figure
pcolor(X(2:nX+1,2:nY+1),Y(2:nX+1,2:nY+1),abs(divV(2:nX+1,2:nY+1)'-divVanal(2:nX+1,2:nY+1)'))
title('divergence error');

%% Compute the potential field for divergence


divPot = zeros(nX+2,nY+2);

divPot = solvePoisson(divPot, divV,xbasis,ybasis,nX,nY);

% errorThresh = 1e-12;
% sorCoeff = 0.9;
% iterMax = 2000;
% % Iteration Loop
% for iter=1:iterMax
%     
% %     fprintf('iteration %d\n',iter);
%     errorMax = -1;
%     
%     % X loop
%     for i=2:nX+1
%         
%         % Y loop
%         for j=2:nY+1
%             % http://ocw.usu.edu/civil_and_environmental_engineering/numerical_methods_in_civil_engineering/LaplaceDirichletTemperature.pdf
%             % Compute hx^2, hy^2, and beta^2, and the denominator for
%             % the old terms
%             
%             alphaX = 2*xbasis(i) - xbasis(i-1) - xbasis(i+1);
%             alphaY = 2*ybasis(j) - ybasis(j-1) - ybasis(j+1);
%             betaX = xbasis(i+1) - xbasis(i-1);
%             betaY = ybasis(j+1) - ybasis(j-1);
%             
%             rX = alphaX/betaX;
%             rY = alphaY/betaY;
%             
%             gammaXSqr = (xbasis(i)-xbasis(i-1))^2 + (xbasis(i+1)-xbasis(i))^2;
%             gammaYSqr = (ybasis(j)-ybasis(j-1))^2 + (ybasis(j+1)-ybasis(j))^2;
%             
%             
%             
%             % Compute the laplacian terms for both components            
%             
%             % Compute the Laplacian discretization
%             potLapla = 0.0;
%             potLapla = potLapla  + gammaYSqr*(divPot(i-1,j)+divPot(i+1,j)+rX*(divPot(i+1,j)-divPot(i-1,j)));
%             potLapla = potLapla  + gammaXSqr*(divPot(i,j-1)+divPot(i,j+1)+rY*(divPot(i,j+1)-divPot(i,j-1)));
%                         
%             
%             % Compute the prefactor on the IJ term
%             % This will divide everything else
%             denom = 2.0*(gammaXSqr + gammaYSqr);
%             
%             % Compute the 
%             
%             
%             % Get the old values
%             divPotOld = divPot(i,j);
%             
%             % Get the new values for the functions velx and vely
%             divPotNewTmp = (potLapla + gammaXSqr*gammaYSqr*divV(i,j))/denom;
%             divPotNew = (1-sorCoeff)*divPotOld + sorCoeff*divPotNewTmp;
%             
%             err = abs(divPotNew-divPotOld)/abs(divPotNew);
%             
%             if(err>errorMax && ~isnan(err))
%                 errorMax = err;
%             end
%             
%             divPot(i,j) = divPotNew;
%             
%         end
%     end
%     
%     
%     % Set boundary conditions
%     divPot(:,end) = divPot(:,end-1); % Neumann on top 
%     divPot(:,1) = divPot(:,2); % Neumann on bottom 
% %     divPot(1,:) = divPot(2,:); % Neumann on left 
% %     divPot(end,:) = divPot(end-1,:); % Neumann on right 
% %     divPot(:,end) = 0;
% %     divPot(:,1) = 0;
%     divPot(1,:) = 0;
%     divPot(end,:) = 0;
%     
%     
%     
% %     divPot(:,end) = divPot(:,end-1); % Neumann on top 
% %     divPot(:,1) = divPot(:,2); % Neumann on bottom 
% %     divPot(1,:) = divPot(2,:) + (divPot(3,:)-divPot(2,:))/(xbasis(3)-xbasis(2))*(xbasis(1)-xbasis(2));
% %     divPot(end,:) = divPot(end-2,:) + (divPot(end-1,:)-divPot(end-2,:))/(xbasis(end-1)-xbasis(end-1))*(xbasis(end)-xbasis(end-2));
%     
%     
%     
% %     divPot(:,end) = divPot(:,2);
% %     divPot(:,1) = divPot(:,end-1);
% %     divPot(1,:) = divPot(end-1,:);
% %     divPot(end,:) = divPot(2,:);
%     
%     
%     
%     if(rem(iter,100)==0 || iter == 1)
%         fprintf('iteration %d\n',iter);
%         fprintf('errorMax: %3.3E\n',errorMax);
%     end
%     
%     if(errorMax<errorThresh)
%         fprintf('Converged! %d %E (%E)\n',iter, errorMax, errorThresh);
%         break;
%     end
%     
% end

%% Remove the gradient of the divergence potential from the velocity components
% X loop
for i=2:nX+1

    % Y loop
    for j=2:nY+1
        
        % http://ocw.usu.edu/civil_and_environmental_engineering/numerical_methods_in_civil_engineering/LaplaceDirichletTemperature.pdf
        
        % Compute the gradients
        betaX = xbasis(i+1) - xbasis(i-1);
        betaY = ybasis(j+1) - ybasis(j-1);
        
        % Compute the source terms
        gradX = (divPot(i+1,j) - divPot(i-1,j))/betaX;
        gradY = (divPot(i,j+1) - divPot(i,j-1))/betaY;
        
        velx(i,j) = velx(i,j) - gradX;
        vely(i,j) = vely(i,j) - gradY;
%         velx(i,j) = velx(i,j) + gradX;
%         vely(i,j) = vely(i,j) + gradY;

    end
end

%% Compute the divergence of the velocity field again

divV(:,:) = 0;

% X loop
for i=2:nX+1

    % Y loop
    for j=2:nY+1
        
        % http://ocw.usu.edu/civil_and_environmental_engineering/numerical_methods_in_civil_engineering/LaplaceDirichletTemperature.pdf
        
        % Compute the gradients
        betaX = xbasis(i+1) - xbasis(i-1);
        betaY = ybasis(j+1) - ybasis(j-1);
        
        % Compute the source terms
        gradX =  (velx(i+1,j) - velx(i-1,j))/betaX;
        gradY =  (vely(i,j+1) - vely(i,j-1))/betaY;
        
        divV(i,j) = gradX + gradY;

    end
end

min(divV(:))
max(divV(:))

%%

figure
hold on;

quiver(X,Y,velx',vely',1)
contour(X,Y,log10(abs(divV')))

%%

figure
potAnal = (1-cos(2*pi*X'))/(2*pi);
pcolor(X,Y,abs(divPot'- potAnal'))
