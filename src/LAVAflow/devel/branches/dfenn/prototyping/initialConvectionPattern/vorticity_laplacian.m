clear all; close all; clc;

nX = 100;
nY = 100;

velx = zeros(nX+2,nY+2);
vely = zeros(nX+2,nY+2);
vorticity = zeros(nX+2,nY+2);



% vorticity(20:40,20:40) = 1.0;





xbnds = [0,1];
ybnds = [0,1];
xbasis = linspace(xbnds(1),xbnds(2),nX+2);
ybasis = linspace(ybnds(1),ybnds(2),nY+2);

hVector = [xbasis(2)-xbasis(1);
           ybasis(2)-ybasis(1)];

[X,Y] = meshgrid(xbasis,ybasis);
Rsqr = (X'-0.5).^2 + (Y'-.5).^2;
       
% vorticity(:,:) = 1000*sin(2*pi*X).*sin(2*pi*Y).*exp(-100*Rsqr);
% vorticity(50:60,50:60) = 1000;

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
        sigma = -1000;
        vorticity = vorticity + amplitude*exp(sigma*Rsqr);    

    end
end

% Circle of constant vorticity
if(1)
    Rsqr = (X'-0.5).^2 + (Y'-0.5).^2;
    la = Rsqr <= 0.01;
    vorticity(la) = 1.0;
    
    Rsqr = (X'-0.75).^2 + (Y'-0.75).^2;
    la = Rsqr <= 0.01;
    vorticity(la) = -1.0;
end

% Compute the gradient of the vorticity field for the source term 
[gradVortX, gradVortY] = gradient2(vorticity, xbasis, ybasis, nX,nY);

       
% Solve the poison equation for the x and y velocities (with dirichlet BC = 0)
% The source term for velx is gradVortX
% The source term for vely is -gradVortY
errorThresh = 1e-5;
iterMax  = 2000;
fprintf('Solving for velx...\n');
velx = solvePoisson(velx,  gradVortX,xbasis,ybasis,nX,nY,'dirichlet',errorThresh,iterMax);
fprintf('Solving for vely...\n');
vely = solvePoisson(vely, -gradVortY,xbasis,ybasis,nX,nY,'dirichlet',errorThresh,iterMax);


% iterMax = 1000;
% 
% errorMax = 1e99;
% errorThresh = 1e-9;
% 
% sorCoeff = 0.9;
% 
% isForward = true;
% % Iteration Loop
% for iter=1:iterMax
%     
% %     fprintf('iteration %d\n',iter);
%     errorMaxX = -1;
%     errorMaxY = -1;
%     
%     isForward = rem(iter,2)==0;
%     
%     offset = 1;
%     if(~isForward)
%         offset = -1;
%     end
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
%             velxLapla = 0.0;
%             velxLapla = velxLapla  + gammaYSqr*(velx(i-1,j)+velx(i+1,j)+rX*(velx(i+1,j)-velx(i-1,j)));
%             velxLapla = velxLapla  + gammaXSqr*(velx(i,j-1)+velx(i,j+1)+rY*(velx(i,j+1)-velx(i,j-1)));
%             
%             
%             velyLapla = 0.0;
%             velyLapla = velyLapla  + gammaYSqr*(vely(i-1,j)+vely(i+1,j)+rX*(vely(i+1,j)-vely(i-1,j)));
%             velyLapla = velyLapla  + gammaXSqr*(vely(i,j-1)+vely(i,j+1)+rY*(vely(i,j+1)-vely(i,j-1)));
%             
%             
%             % Compute the prefactor on the IJ term
%             % This will divide everything else
%             denom = 2.0*(gammaXSqr + gammaYSqr);
%             
%             % Compute the source terms
%             velxSource =  (vorticity(i,j+1) - vorticity(i,j-1))/betaX;
%             velySource = -(vorticity(i+1,j) - vorticity(i-1,j))/betaY;
%             
%             
%             % Get the old values
%             velxOld = velx(i,j);
%             velyOld = vely(i,j);
%             
%             % Get the new values for the functions velx and vely
%             velxNewTmp = (velxLapla + gammaXSqr*gammaYSqr*velxSource)/denom;
%             velyNewTmp = (velyLapla + gammaXSqr*gammaYSqr*velySource)/denom;
%             
%             velxNew = (1-sorCoeff)*velxOld + sorCoeff*velxNewTmp;
%             velyNew = (1-sorCoeff)*velyOld + sorCoeff*velyNewTmp;
%             
%             errorX = abs(velxNew-velxOld)/abs(velxNew);
%             errorY = abs(velyNew-velyOld)/abs(velyNew);
%             
% %             errorX = abs(velxNew-velxOld);
% %             errorY = abs(velyNew-velyOld);
%             
%             if(errorX>errorMaxX)
%                 errorMaxX = errorX;
%             end
%             if(errorY>errorMaxY)
%                 errorMaxY = errorY;
%             end
%             
%             velx(i,j) = velxNew;
%             vely(i,j) = velyNew;
%             
% %             
% %             % Compute the contribution from the source terms
% %             sourceTerms = gammaXSqr*gammaYSqr*source(i,j);
% %             
% %             % Compute the denominator
% %             denom = 2.0*(gammaXSqr + gammaYSqr);
% %             
% %             
% %             potIJnew = (dens(i,j)*derivTerms + sourceTerms - gradTerms)/(dens(i,j)*denom);
% %             
% %             potIJnewTmp = (1-sorCoeff)*pot(i,j) + sorCoeff*potIJnew;
% %             
% %             error = abs((pot(i,j)-potIJnewTmp));
% %             if(error>errorMax)
% %                 errorMax = error;
% %             end
% %             
% %             pot(i,j) = potIJnewTmp;
%             
%         end
%     end
%     
%     
%     % Set boundary conditions
% %     velx(:,end) = velx(:,end-1); % Neumann on top 
% %     velx(:,1) = velx(:,2); % Dirichlet on bottom 
% %     velx(1,:) = velx(2,:); % Neumann on left 
% %     velx(end,:) = velx(end-1,:); % Neumann on right 
% % %     velx(end,:) = 0; % Neumann on right 
% %     
% %     
% %     vely(:,end) = vely(:,end-1); % Neumann on top 
% %     vely(:,1) = 0; % Dirichlet on bottom 
% %     vely(1,:) = vely(2,:); % Neumann on left 
% %     vely(end,:) = vely(end-1,:); % Neumann on right 
%     
%     %====================================
%     
% %     velx(:,end) = 0; % Neumann on top 
% %     velx(:,1) = velx(:,2); % Dirichlet on bottom 
% %     velx(1,:) = 0; % Neumann on left 
% %     velx(end,:) = 0; % Neumann on right 
% %     
% %     
% %     vely(:,end) = vely(:,end-1); % Neumann on top 
% %     vely(:,1) = 0; % Dirichlet on bottom 
% %     vely(1,:) = vely(2,:); % Neumann on left 
% %     vely(end,:) = vely(end-1,:); % Neumann on right 
%     
%     %====================================
%     
%     velx(:,end) = 0;
%     velx(:,1) = 0; % Dirichlet on bottom 
%     velx(1,:) = 0;
%     velx(end,:) = 0;
%     
%     vely(:,end) = 0;
%     vely(:,1) = 0;
%     vely(1,:) = 0;
%     vely(end,:) = 0;
%     
%     
%     
%     velx(:,end) = velx(:,2);
%     velx(:,1) = velx(:,end-1); 
%     velx(1,:) = velx(end-1,:);
%     velx(end,:) = velx(2,:);
%     vely(:,end) = vely(:,2);
%     vely(:,1) = vely(:,end-1); 
%     vely(1,:) = vely(end-1,:);
%     vely(end,:) = vely(2,:);
%     
%     if(rem(iter,100)==0 || iter == 1)
%         fprintf('iteration %d\n',iter);
%         fprintf('errorMaxX: %3.3E\n',errorMaxX);
%         fprintf('errorMaxY: %3.3E\n',errorMaxY);
%     end
%     
%     if(errorMaxX<errorThresh && errorMaxY<errorThresh)
%         fprintf('Converged!\n');
%         break;
%     end
%     
% end
% 
% 
% figure
% hold on;
% 
% % [gradX,gradY] = gradient(pot,hVector(1),hVector(2));
% [X,Y] = meshgrid(xbasis,ybasis);
% 
% [curlZ,~] = curl(X,Y,velx',vely');
% 
% V = sqrt(velx.^2+vely.^2);
% la = V>=0.1*max(V(:));
% 
% contour(X,Y,vorticity')
% quiver(X,Y,velx',vely',3)
% 
% axis equal

%% Compute the divergence of the velocity field

divV = divergence2(velx,vely,xbasis,ybasis,nY,nX);

fprintf('Initial divergence:\n');
fprintf('\tMin: %E\n',min(divV(:)));
fprintf('\tMax: %E\n',max(divV(:)));

% divV = zeros(nX+2,nY+2);
% 
% % X loop
% for i=2:nX+1
% 
%     % Y loop
%     for j=2:nY+1
%         
%         % http://ocw.usu.edu/civil_and_environmental_engineering/numerical_methods_in_civil_engineering/LaplaceDirichletTemperature.pdf
%         
%         % Compute the gradients
%         betaX = xbasis(i+1) - xbasis(i-1);
%         betaY = ybasis(j+1) - ybasis(j-1);
%         
%         % Compute the source terms
%         gradX =  (velx(i+1,j) - velx(i-1,j))/betaX;
%         gradY =  (vely(i,j+1) - vely(i,j-1))/betaY;
%         
%         divV(i,j) = gradX + gradY;
% 
%     end
% end
% 
% % Compute the divergence on the boundary via forward/backward differences
% 
% % Left boundary
% for i=1:1
%     for j=2:nY+1
%         
%         % Compute the gradients
%         betaX = xbasis(i+2) - xbasis(i);
%         betaY = ybasis(j+1) - ybasis(j-1);
%         
%         % Compute the source terms
%         gradX =  (-velx(i+2,j) + 4.0*velx(i+1,j) - 3.0*velx(i,j))/betaX;
%         gradY =  (vely(i,j+1) - vely(i,j-1))/betaY;
%         divV(i,j) = gradX + gradY;
%     end
% end
% 
% 
% % Right boundary
% for i=nX+2:nX+2
%     for j=2:nY+1
%         
%         % Compute the gradients
%         betaX = xbasis(i) - xbasis(i-2);
%         betaY = ybasis(j+1) - ybasis(j-1);
%         
%         % Compute the source terms
%         gradX =  (velx(i-2,j) - 4.0*velx(i-1,j) + 3.0*velx(i,j))/betaX;
%         gradY =  (vely(i,j+1) - vely(i,j-1))/betaY;
%         divV(i,j) = gradX + gradY;
%     end
% end
% 
% 
% % Bottom boundary
% for i=2:nX+1
%     for j=1:1
%         
%         % Compute the gradients
%         betaX = xbasis(i+1) - xbasis(i-1);
%         betaY = ybasis(j+2) - ybasis(j);
%         
%         % Compute the source terms
%         gradX =  (velx(i+1,j) - velx(i-1,j))/betaX;
%         gradY =  (-vely(i,j+2) + 4.0*vely(i,j+1) - 3.0*vely(i,j))/betaY;
%         divV(i,j) = gradX + gradY;
%     end
% end
% 
% 
% % Top boundary
% for i=2:nX+1
%     for j=nY+2:nY+2
%         
%         % Compute the gradients
%         betaX = xbasis(i+1) - xbasis(i-1);
%         betaY = ybasis(j) - ybasis(j-2);
%         
%         % Compute the source terms
%         gradX =  (velx(i+1,j) - velx(i-1,j))/betaX;
%         gradY =  (vely(i,j-2) - 4.0*vely(i,j-1) + 3.0*vely(i,j))/betaY;
%         divV(i,j) = gradX + gradY;
%     end
% end
% 
% min(divV(:))
% max(divV(:))
% 
% V = sqrt(velx.^2+vely.^2);
% divVImportance = hVector(1)*divV./max(V(:));
% min(abs(divVImportance(divVImportance~=Inf)))
% max(abs(divVImportance(divVImportance~=Inf)))
% max(V(:))

%% Compute the potential field for divergence


fprintf('Solving for the potential field...\n');
errorThresh = 1e-6;
iterMax = 20000;
divPot = zeros(nX+2,nY+2);
% divPot = solvePoisson(divPot, divV,xbasis,ybasis,nX,nY,'periodic',errorThresh,iterMax);
% divPot = solvePoisson(divPot, divV,xbasis,ybasis,nX,nY,'dirichlet',errorThresh,iterMax);
divPot = solvePoisson(divPot, divV,xbasis,ybasis,nX,nY,'outflow',errorThresh,iterMax);


% divPot = zeros(nX+2,nY+2);
% errorThresh = 1e-12;
% sorCoeff = 0.5;
% iterMax = 1000;
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
%             err = abs(divPotNew-divPotOld);
%             
%             if(err>errorMax)
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
%     divPot(1,:) = divPot(2,:); % Neumann on left 
%     divPot(end,:) = divPot(end-1,:); % Neumann on right 
% %     divPot(:,end) = 0;
% %     divPot(:,1) = 0;
% %     divPot(1,:) = 0;
% %     divPot(end,:) = 0;
%     
% 
%     divPot(:,end) = divPot(:,2);
%     divPot(:,1) = divPot(:,end-1);
%     divPot(1,:) = divPot(end-1,:);
%     divPot(end,:) = divPot(2,:);
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


fprintf('Removing compressive component...\n');
[gradPotX, gradPotY] = gradient2(divPot, xbasis, ybasis, nX,nY);

velx = velx - gradPotX;
vely = vely - gradPotY;


% % X loop
% for i=2:nX+1
% 
%     % Y loop
%     for j=2:nY+1
%         
%         % http://ocw.usu.edu/civil_and_environmental_engineering/numerical_methods_in_civil_engineering/LaplaceDirichletTemperature.pdf
%         
%         % Compute the gradients
%         betaX = xbasis(i+1) - xbasis(i-1);
%         betaY = ybasis(j+1) - ybasis(j-1);
%         
%         % Compute the source terms
%         gradX = (divPot(i+1,j) - divPot(i-1,j))/betaX;
%         gradY = (divPot(i,j+1) - divPot(i,j-1))/betaY;
%         
%         velx(i,j) = velx(i,j) - gradX;
%         vely(i,j) = vely(i,j) - gradY;
% %         velx(i,j) = velx(i,j) + gradX;
% %         vely(i,j) = vely(i,j) + gradY;
% 
%     end
% end



%% Compute the divergence of the velocity field again


divV = divergence2(velx,vely,xbasis,ybasis,nY,nX);

fprintf('Final divergence:\n');
fprintf('\tMin: %E\n',min(divV(:)));
fprintf('\tMax: %E\n',max(divV(:)));

% divV(:,:) = 0;
% 
% % X loop
% for i=2:nX+1
% 
%     % Y loop
%     for j=2:nY+1
%         
%         % http://ocw.usu.edu/civil_and_environmental_engineering/numerical_methods_in_civil_engineering/LaplaceDirichletTemperature.pdf
%         
%         % Compute the gradients
%         betaX = xbasis(i+1) - xbasis(i-1);
%         betaY = ybasis(j+1) - ybasis(j-1);
%         
%         % Compute the source terms
%         gradX =  (velx(i+1,j) - velx(i-1,j))/betaX;
%         gradY =  (vely(i,j+1) - vely(i,j-1))/betaY;
%         
%         divV(i,j) = gradX + gradY;
% 
%     end
% end
% 
% % Compute the divergence on the boundary via forward/backward differences
% 
% % Left boundary
% for i=1:1
%     for j=2:nY+1
%         
%         % Compute the gradients
%         betaX = xbasis(i+2) - xbasis(i);
%         betaY = ybasis(j+1) - ybasis(j-1);
%         
%         % Compute the source terms
%         gradX =  (-velx(i+2,j) + 4.0*velx(i+1,j) - 3.0*velx(i,j))/betaX;
%         gradY =  (vely(i,j+1) - vely(i,j-1))/betaY;
%         divV(i,j) = gradX + gradY;
%     end
% end
% 
% 
% % Right boundary
% for i=nX+2:nX+2
%     for j=2:nY+1
%         
%         % Compute the gradients
%         betaX = xbasis(i) - xbasis(i-2);
%         betaY = ybasis(j+1) - ybasis(j-1);
%         
%         % Compute the source terms
%         gradX =  (velx(i-2,j) - 4.0*velx(i-1,j) + 3.0*velx(i,j))/betaX;
%         gradY =  (vely(i,j+1) - vely(i,j-1))/betaY;
%         divV(i,j) = gradX + gradY;
%     end
% end
% 
% 
% % Bottom boundary
% for i=2:nX+1
%     for j=1:1
%         
%         % Compute the gradients
%         betaX = xbasis(i+1) - xbasis(i-1);
%         betaY = ybasis(j+2) - ybasis(j);
%         
%         % Compute the source terms
%         gradX =  (velx(i+1,j) - velx(i-1,j))/betaX;
%         gradY =  (-vely(i,j+2) + 4.0*vely(i,j+1) - 3.0*vely(i,j))/betaY;
%         divV(i,j) = gradX + gradY;
%     end
% end
% 
% 
% % Top boundary
% for i=2:nX+1
%     for j=nY+2:nY+2
%         
%         % Compute the gradients
%         betaX = xbasis(i+1) - xbasis(i-1);
%         betaY = ybasis(j) - ybasis(j-2);
%         
%         % Compute the source terms
%         gradX =  (velx(i+1,j) - velx(i-1,j))/betaX;
%         gradY =  (vely(i,j-2) - 4.0*vely(i,j-1) + 3.0*vely(i,j))/betaY;
%         divV(i,j) = gradX + gradY;
%     end
% end
% 
% min(divV(:))
% max(divV(:))
            
%%
figure
hold on;

% [gradX,gradY] = gradient(pot,hVector(1),hVector(2));
[X,Y] = meshgrid(xbasis,ybasis);

V = sqrt(velx.^2+vely.^2);
la = V>=0.1*max(V(:));

contour(X,Y,vorticity')
quiver(X,Y,velx',vely',3)

axis equal


%%



figure
h=pcolor(X,Y,log10(abs(divV')));
set(h,'edgealpha',0);

            

            
            


