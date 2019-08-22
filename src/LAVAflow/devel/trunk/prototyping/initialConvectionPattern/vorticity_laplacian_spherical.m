clear all; close all; clc;

nX = 100;
nY = 100;

velr = zeros(nX+2,nY+2);
velp = zeros(nX+2,nY+2);
vorticity = zeros(nX+2,nY+2);



% vorticity(20:40,20:40) = 1.0;


% =============================================
%              Domain information
% =============================================

xbnds = [1e6,1e8];
ybnds = [pi/4,3*pi/4];
xbasis = linspace(xbnds(1),xbnds(2),nX+2);
ybasis = linspace(ybnds(1),ybnds(2),nY+2);

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

        amplitude = (-1)^i*1e8;
        sigma = -100;
        
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
dens = ones(nX+2,nY+2);
gradDensR = zeros(nX+2,nY+2);
gradDensP = zeros(nX+2,nY+2);

if(0)
    dens(:,:) = 1.0;
end

if(0)
    Lr = 1;
    Lp = 1;
    Amp = 1;
    dens = Amp*(sin(2*pi*(PHI'-ybnds(1))/(ybnds(2)-ybnds(1))*Lp).*sin(2*pi*R'/Lr)+1.01)/2;
end


if(1)
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

% Interior
for i=2:nX+1
    r = xbasis(i);
    for j=2:nY+1
        
        phi = ybasis(j);
        betaX = xbasis(i+1) - xbasis(i-1);
        betaY = ybasis(j+1) - ybasis(j-1);
        
        gradDensR(i,j) = (dens(i+1,j)-dens(i-1,j))/betaX;
        gradDensP(i,j) = (1.0/r)*(dens(i,j+1)-dens(i,j-1))/betaY;
        
    end
end

% Bottom boundary
for i=1:1
    r = xbasis(i);
    for j=2:nY+1
        
        betaX = xbasis(i+1) - xbasis(i);
        betaY = ybasis(j+1) - ybasis(j-1);
        gradDensR(i,j) = (dens(i+1,j)-dens(i,j))/betaX;
        gradDensP(i,j) = (1.0/r)*(dens(i,j+1)-dens(i,j-1))/betaY;
        
    end
end


% Top boundary
for i=nX+2:nX+2
    r = xbasis(i);
    for j=2:nY+1
        
        betaX = xbasis(i) - xbasis(i-1);
        betaY = ybasis(j+1) - ybasis(j-1);
        gradDensR(i,j) = (dens(i,j)-dens(i-1,j))/betaX;
        gradDensP(i,j) = (1.0/r)*(dens(i,j+1)-dens(i,j-1))/betaY;
        
    end
end


% Left boundary
for i=2:nX+1
    r = xbasis(i);
    for j=1:1
        
        betaX = xbasis(i+1) - xbasis(i-1);
        betaY = ybasis(j+1) - ybasis(j);
        gradDensR(i,j) = (dens(i+1,j)-dens(i-1,j))/betaX;
        gradDensP(i,j) = (1.0/r)*(dens(i,j+1)-dens(i,j))/betaY;
        
    end
end


% Right boundary
for i=2:nX+1
    r = xbasis(i);
    for j=nY+2:nY+2
        
        betaX = xbasis(i+1) - xbasis(i-1);
        betaY = ybasis(j) - ybasis(j-1);
        gradDensR(i,j) = (dens(i+1,j)-dens(i-1,j))/betaX;
        gradDensP(i,j) = (1.0/r)*(dens(i,j)-dens(i,j-1))/betaY;
        
    end
end
        

% % Right boundary
% gradDensR(end,2:nY+1) = (dens(end,3:nY+2)-dens(end,1:nY))./(xbasis(3:nY+2)-xbasis(1:nY));
% gradDensP(end,2:nY+1) = (dens(end,2:nY+1)-dens(end-1,2:nY+1))./(ybasis(end)-ybasis(end-1));
% 
% 
% % Left boundary
% gradDensR(1,2:nY+1) = (dens(1,3:nY+2)-dens(1,1:nY))./(xbasis(3:nY+2)-xbasis(1:nY));
% gradDensP(1,2:nY+1) = (dens(2,2:nY+1)-dens(1,2:nY+1))./(ybasis(2)-ybasis(1));
% 
% 
% % Upper boundary
% gradDensR(2:nX+1,end) = (dens(2:nX+1,end)-dens(2:nX+1,end-1))./(xbasis(end)-xbasis(end-1));
% gradDensP(2:nX+1,end) = (dens(3:nX+2,end)-dens(1:nX,end))./(ybasis(3:nY+2)-ybasis(1:nY));
% 
% 
% % Lower boundary
% gradDensR(2:nX+1,1) = (dens(2:nX+1,2)-dens(2:nX+1,1))./(xbasis(2)-xbasis(1));
% gradDensP(2:nX+1,1) = (dens(3:nX+2,1)-dens(1:nX,1))./(ybasis(3:nY+2)-ybasis(1:nY));



% =============================================
%           Solve for velocity field
% =============================================
divergenceSource = zeros(nX+2,nY+2);

iterMax = 1000;

errorMax = 1e99;
errorThresh = 1e-9;

sorCoeff = 0.5;

isForward = true;
% Iteration Loop
for iter=1:iterMax
    
%     fprintf('iteration %d\n',iter);
    errorMaxX = -1;
    errorMaxY = -1;
    
    isForward = rem(iter,2)==0;
    
    offset = 1;
    if(~isForward)
        offset = -1;
    end
    
    % Compute the term dot(u,gradDens)/dens
    divergenceSource = (velr.*gradDensR + velp.*gradDensP)./dens;
    
    % X loop
    for i=2:nX+1
        
        r = xbasis(i);
        
        % Y loop
        for j=2:nY+1
            
            phi = ybasis(j);
            
            % http://ocw.usu.edu/civil_and_environmental_engineering/numerical_methods_in_civil_engineering/LaplaceDirichletTemperature.pdf
            % Compute hx^2, hy^2, and beta^2, and the denominator for
            % the old terms
            
            alphaX = 2*xbasis(i) - xbasis(i-1) - xbasis(i+1);
            alphaY = 2*ybasis(j) - ybasis(j-1) - ybasis(j+1);
            betaX = xbasis(i+1) - xbasis(i-1);
            betaY = ybasis(j+1) - ybasis(j-1);
            
            rX = alphaX/betaX;
            rY = alphaY/betaY;
            
            gammaXSqr = (xbasis(i)-xbasis(i-1))^2 + (xbasis(i+1)-xbasis(i))^2;
            gammaYSqr = (ybasis(j)-ybasis(j-1))^2 + (ybasis(j+1)-ybasis(j))^2;
            
            
            
            % Compute the laplacian terms for both components            
            
                        
            % Compute the Laplacian discretization
            velrLapla = 0.0;
            velrLapla = velrLapla  + gammaYSqr*(velr(i-1,j)+velr(i+1,j)+rX*(velr(i+1,j)-velr(i-1,j)));
            velrLapla = velrLapla  + gammaXSqr*(velr(i,j-1)+velr(i,j+1)+rY*(velr(i,j+1)-velr(i,j-1)))/xbasis(i)^2;
            velrLapla = velrLapla  + 2.0*gammaXSqr*gammaYSqr*(velr(i+1,j)-velr(i-1,j))/(betaX*xbasis(i));
            velrLapla = velrLapla  + gammaXSqr*gammaYSqr*(velr(i,j+1)-velr(i,j-1))/(betaY*xbasis(i)^2)*cot(ybasis(j));
            
            
            
            
            velpLapla = 0.0;
            velpLapla = velpLapla  + gammaYSqr*(velp(i-1,j)+velp(i+1,j)+rX*(velp(i+1,j)-velp(i-1,j)));
            velpLapla = velpLapla  + gammaXSqr*(velp(i,j-1)+velp(i,j+1)+rY*(velp(i,j+1)-velp(i,j-1)))/xbasis(i)^2;
            velpLapla = velpLapla  + 2.0*gammaXSqr*gammaYSqr*(velp(i+1,j)-velp(i-1,j))/(betaX*xbasis(i));
            velpLapla = velpLapla  + gammaXSqr*gammaYSqr*(velp(i,j+1)-velp(i,j-1))/(betaY*xbasis(i)^2)*cot(ybasis(j));
            
            
            % Compute the source terms for the velr Poisson equation
            velrSource = 0.0;
            velrSource = velrSource - (divergenceSource(i+1,j)-divergenceSource(i-1,j))/betaX;
            velrSource = velrSource + (1.0/r)*(vorticity(i,j+1) - vorticity(i,j-1))/betaY;
            velrSource = velrSource + 2.0*(velp(i,j+1)-velp(i,j-1))/betaY/r^2;
            velrSource = velrSource + cot(phi)*velp(i,j)/r^2;
            velrSource = velrSource - cot(phi)*( (velp(i+1,j)-velp(i-1,j))/betaX - (velr(i,j+1)-velr(i,j-1))/betaY/r );
            
            velpSource = 0.0;
            velpSource = velpSource - (1.0/r)*(divergenceSource(i,j+1)-divergenceSource(i,j-1))/betaY;
            velpSource = velpSource - (vorticity(i+1,j)-vorticity(i-1,j))/betaX;
            velpSource = velpSource + (1.0/r)*(velp(i+1,j)-velp(i-1,j))/betaX;
            velpSource = velpSource + (1.0/r^2)*(1+csc(phi)^2)*velp(i,j);
            velpSource = velpSource + (3.0/r^2)*(velr(i+1,j)-velr(i-1,j))/betaX;
            
            % Compute the prefactor on the IJ term
            % This will divide everything else
            denom = 2.0*(gammaXSqr + gammaYSqr);
                      
            
            % Get the old values
            velrOld = velr(i,j);
            velpOld = velp(i,j);
            
            % Get the new values for the functions velr and velp
            velrNewTmp = (velrLapla + gammaXSqr*gammaYSqr*velrSource)/denom;
            velpNewTmp = (velpLapla + gammaXSqr*gammaYSqr*velpSource)/denom;
            
            velrNew = (1-sorCoeff)*velrOld + sorCoeff*velrNewTmp;
            velpNew = (1-sorCoeff)*velpOld + sorCoeff*velpNewTmp;
            
            errorX = abs(velrNew-velrOld)/abs(velrNew);
            errorY = abs(velpNew-velpOld)/abs(velpNew);
            
            if(errorX>errorMaxX)
                errorMaxX = errorX;
            end
            if(errorY>errorMaxY)
                errorMaxY = errorY;
            end
            
            velr(i,j) = velrNew;
            velp(i,j) = velpNew;
                        
        end
    end
    
    
    % Set boundary conditions
%     velr(:,end) = velr(:,end-1); % Neumann on top 
%     velr(:,1) = velr(:,2); % Dirichlet on bottom 
%     velr(1,:) = velr(2,:); % Neumann on left 
%     velr(end,:) = velr(end-1,:); % Neumann on right 
% %     velr(end,:) = 0; % Neumann on right 
%     
%     
%     velp(:,end) = velp(:,end-1); % Neumann on top 
%     velp(:,1) = 0; % Dirichlet on bottom 
%     velp(1,:) = velp(2,:); % Neumann on left 
%     velp(end,:) = velp(end-1,:); % Neumann on right 
    
    %====================================
    
%     velr(:,end) = 0; % Neumann on top 
%     velr(:,1) = velr(:,2); % Dirichlet on bottom 
%     velr(1,:) = 0; % Neumann on left 
%     velr(end,:) = 0; % Neumann on right 
%     
%     
%     velp(:,end) = velp(:,end-1); % Neumann on top 
%     velp(:,1) = 0; % Dirichlet on bottom 
%     velp(1,:) = velp(2,:); % Neumann on left 
%     velp(end,:) = velp(end-1,:); % Neumann on right 
    
    %====================================
    
    velr(:,end) = 0;
    velr(:,1) = 0; % Dirichlet on bottom 
    velr(1,:) = 0;
    velr(end,:) = 0;
    
    velp(:,end) = 0;
    velp(:,1) = 0;
    velp(1,:) = 0;
    velp(end,:) = 0;
    
    if(rem(iter,10)==0 || iter == 1)
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
        
        
        gradX = (dens(i+1,j)*velr(i+1,j)-dens(i-1,j)*velr(i-1,j))/betaX;
        gradY = (dens(i,j+1)*velp(i,j+1)-dens(i,j-1)*velp(i,j-1))/(betaY*xbasis(i));
        
        divV(i,j) = gradX + gradY + 2.0*dens(i,j)*velr(i,j)/xbasis(i) + cot(ybasis(j))*dens(i,j)*velp(i,j)/xbasis(i);
        

    end
end


%% Compute the potential field for divergence


divPot = zeros(nX+2,nY+2);


errorThresh = 1e-9;
sorCoeff = 0.75;
% Iteration Loop
for iter=1:iterMax
    
%     fprintf('iteration %d\n',iter);
    errorMax = -1;
    
    % X loop
    for i=2:nX+1
        
        % Y loop
        for j=2:nY+1
            % http://ocw.usu.edu/civil_and_environmental_engineering/numerical_methods_in_civil_engineering/LaplaceDirichletTemperature.pdf
            % Compute hx^2, hy^2, and beta^2, and the denominator for
            % the old terms
            
            alphaX = 2*xbasis(i) - xbasis(i-1) - xbasis(i+1);
            alphaY = 2*ybasis(j) - ybasis(j-1) - ybasis(j+1);
            betaX = xbasis(i+1) - xbasis(i-1);
            betaY = ybasis(j+1) - ybasis(j-1);
            
            rX = alphaX/betaX;
            rY = alphaY/betaY;
            
            gammaXSqr = (xbasis(i)-xbasis(i-1))^2 + (xbasis(i+1)-xbasis(i))^2;
            gammaYSqr = (ybasis(j)-ybasis(j-1))^2 + (ybasis(j+1)-ybasis(j))^2;
            
            
            
            % Compute the laplacian terms for both components            
            
            % Compute the Laplacian discretization
            potLapla = 0.0;
            potLapla = potLapla  + gammaYSqr*(divPot(i-1,j)+divPot(i+1,j)+rX*(divPot(i+1,j)-divPot(i-1,j)));
            potLapla = potLapla  + gammaXSqr*(divPot(i,j-1)+divPot(i,j+1)+rY*(divPot(i,j+1)-divPot(i,j-1)))/xbasis(i)^2;
            potLapla = potLapla  + 2.0*gammaXSqr*gammaYSqr*(divPot(i+1,j)-divPot(i-1,j))/(betaX*xbasis(i));
            potLapla = potLapla  + gammaXSqr*gammaYSqr*(divPot(i,j+1)-divPot(i,j-1))/(betaY*xbasis(i)^2)*cot(ybasis(j));
                        
            
            % Compute the prefactor on the IJ term
            % This will divide everything else
            denom = 2.0*(gammaXSqr + gammaYSqr);
            
            % Compute the 
            
            
            % Get the old values
            divPotOld = divPot(i,j);
            
            % Get the new values for the functions velx and vely
            divPotNewTmp = (potLapla + gammaXSqr*gammaYSqr*divV(i,j))/denom;
            divPotNew = (1-sorCoeff)*divPotOld + sorCoeff*divPotNewTmp;
            
            err = abs(divPotNew-divPotOld)/abs(divPotNew);
            
            if(err>errorMax)
                errorMax = err;
            end
            
            divPot(i,j) = divPotNew;
            
        end
    end
    
    
    % Set boundary conditions
    divPot(:,end) = divPot(:,end-1); % Neumann on top 
    divPot(:,1) = divPot(:,2); % Neumann on bottom 
    divPot(1,:) = divPot(2,:); % Neumann on left 
    divPot(end,:) = divPot(end-1,:); % Neumann on right 
    
    if(rem(iter,100)==0 || iter == 1)
        fprintf('iteration %d\n',iter);
        fprintf('errorMax: %3.3E\n',errorMax);
    end
    
    if(errorMax<errorThresh)
        fprintf('Converged! %d %E (%E)\n',iter, errorMax, errorThresh);
        break;
    end
    
end


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
        gradY = ((divPot(i,j+1) - divPot(i,j-1))/betaY)/xbasis(i);
        
        velr(i,j) = velr(i,j) - gradX;
        velp(i,j) = velp(i,j) - gradY;

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
        
        gradX = (dens(i+1,j)*velr(i+1,j)-dens(i-1,j)*velr(i-1,j))/betaX;
        gradY = (dens(i,j+1)*velp(i,j+1)-dens(i,j-1)*velp(i,j-1))/(betaY*xbasis(i));
        divV(i,j) = gradX + gradY + 2.0*dens(i,j)*velr(i,j)/xbasis(i) + cot(ybasis(j))*dens(i,j)*velp(i,j)/xbasis(i);
    end
end

min(divV(:))
max(divV(:))
            
%% Convert vectors into cartesian
velx = zeros(nX+2,nY+2);
vely = zeros(nX+2,nY+2);

% X loop
for i=2:nX+1

    % Y loop
    for j=2:nY+1
        
        betaX = xbasis(i+1) - xbasis(i-1);
        betaY = ybasis(j+1) - ybasis(j-1);
        
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
h=pcolor(X,Y,log10(abs(divPot')));
set(h,'edgealpha',0);

min(divV(:))
max(divV(:))
            

            
            


