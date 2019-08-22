clear all; close all; clc; pause(0.1);

nX = 100;
nY = 100;

pot = zeros(nX+2,nY+2);
dens = ones(nX+2,nY+2);
source = zeros(nX+2,nY+2);


% xbnds = [1e8,1e9];
xbnds = [1,100];
ybnds = [pi/4,3*pi/4];
radExtent = diff(xbnds);
radOffset = radExtent/4;
phiOffset = 0.3;

xbasis = linspace(xbnds(1),xbnds(2),nX+2);
% xbasis = logspace(log10(xbnds(1)),log10(xbnds(2)),nX+2);

ybasis = linspace(ybnds(1),ybnds(2),nY+2);

hVector = [xbasis(2)-xbasis(1);
           ybasis(2)-ybasis(1)];


[R,PHI] = meshgrid(xbasis,ybasis);
X = R.*cos(PHI);
Y = R.*sin(PHI);
Rsqr = (X-0).^2 + (Y-1.5).^2;
       
% Egg crate density
Lr = 2.0;
Lp = 2;
% dens = sin(2*pi*(PHI'-ybnds(1))/(ybnds(2)-ybnds(1))*Lp).*sin(2*pi*R'/Lr)+1.1;
dens(:,:) = 1.0;

       
% Random Spherical Gaussians
if(1)
    xmin = 1e6;
    xmax = 9e7;
    source(:,:) = 0.0;
    nGauss = 1000;
    sigmaArr = zeros(nGauss,1);
    for i=1:nGauss
        posSphere = [xbnds(1)+radOffset, ybnds(1)+phiOffset] + rand(1,2).*[diff(xbnds)-2*radOffset,diff(ybnds)-2*phiOffset];
        pos = [posSphere(1)*cos(posSphere(2)), posSphere(1)*sin(posSphere(2))];
        Rsqr = (X'-pos(1)).^2 + (Y'-pos(2)).^2;

    %     amplitude = (-1)^i*((posSphere(1)-xbnds(1))/(xbnds(2)-xbnds(1)-offset)*1000 + 0);

        amplitude = (-1)^i;
        sigma = -10;
        if(posSphere(1)>1e9)
            amplitude = amplitude*1;
        else
            amplitude = amplitude*radExtent^2;
            sigma = -100;
        end


        sigma = -10000 + (posSphere(1)-xbnds(1))/(xbnds(2)-xbnds(1)-radOffset)*10000;
        source = source + amplitude*exp(sigma*Rsqr/radExtent^2);    

        sigmaArr(i) = sigma;
    end
end

% Cylinders
if(0)
    source(:,:) = 0.0;
    nCyl = 50;
    offset = 0.3;
    for i=1:nCyl
        posSphere = [xbnds(1)+offset, ybnds(1)+offset] + rand(1,2).*[diff(xbnds)-2*offset,diff(ybnds)-2*offset];
        pos = [posSphere(1)*cos(posSphere(2)), posSphere(1)*sin(posSphere(2))];
        Rsqr = (X'-pos(1)).^2 + (Y'-pos(2)).^2;

        if(posSphere(1)<5e7)
            i = i-1;
            continue;
        end
        
        la = Rsqr<5e14;

        source(la) = source(la) + 1;
        
    end
end


source(:,1) = 0.0;
source(:,end) = 0.0;
source(1,:) = 0.0;
source(end,:) = 0.0;

% Initial guess
% pot(:,:) = source(:,:);

%%       
% Dirichlet of pot=1 for all x at the upper boundary
% pot(:,end) = 1.0;

iterMax = 20;

errorMax = 1e99;
errorThresh = 1e-9;

sorCoeff = 0.1;

% Iteration Loop
for iter=1:iterMax
    
    fprintf('iteration %d\n',iter);
    errorMax = 0;
    
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
            
                        
            
            % Compute the Laplacian discretization
            derivTerms = 0.0;
            derivTerms = derivTerms + gammaYSqr*(pot(i-1,j)+pot(i+1,j)+rX*(pot(i+1,j)-pot(i-1,j)));
            derivTerms = derivTerms + gammaXSqr*(pot(i,j-1)+pot(i,j+1)+rY*(pot(i,j+1)-pot(i,j-1)))/xbasis(i)^2;
            derivTerms = derivTerms + 2.0*gammaXSqr*gammaYSqr*(pot(i+1,j)-pot(i-1,j))/(betaX*xbasis(i));
            derivTerms = derivTerms + gammaXSqr*gammaYSqr*(pot(i,j+1)-pot(i,j-1))/(betaY*xbasis(i)^2)*cot(ybasis(j));
            
            % Compute the contribution from the gradDens_dot_gradPot term
            gradDens = zeros(3,1);
            gradDens(1) = (dens(i+1,j) - dens(i-1,j))/betaX;
            gradDens(2) = (dens(i,j+1) - dens(i,j-1))/(betaY*xbasis(i));
            
            gradPot = zeros(3,1);
            gradPot(1) = (pot(i+1,j) - pot(i-1,j))/betaX;
            gradPot(2) = (pot(i,j+1) - pot(i,j-1))/(betaY*xbasis(i));
            
            gradTerms = gammaXSqr*gammaYSqr*dot(gradDens,gradPot);
            
            % Compute the contribution from the source terms
            sourceTerms = gammaXSqr*gammaYSqr*source(i,j);
            
            % Compute the denominator
            denom = 2.0*(gammaXSqr + gammaYSqr);
            
            
            potIJnew = (dens(i,j)*derivTerms + sourceTerms - gradTerms)/(dens(i,j)*denom);
            
            potIJnewTmp = (1-sorCoeff)*pot(i,j) + sorCoeff*potIJnew;
            
            error = abs((pot(i,j)-potIJnewTmp));
            if(error>errorMax)
                errorMax = error;
            end
            
            pot(i,j) = potIJnewTmp;
            
        end
    end
    
    fprintf('errorMax: %3.3E\n',errorMax);
    
    if(errorMax<errorThresh)
        fprintf('Converged!\n');
        break;
    end
    
end
            
            
%% Compute velocity and divergence

velr = zeros(nX+2,nY+2);
velp = zeros(nX+2,nY+2);
velx = zeros(nX+2,nY+2);
vely = zeros(nX+2,nY+2);
divV = zeros(nX+2,nY+2);

% X loop
for i=2:nX+1

    % Y loop
    for j=2:nY+1
        
        betaX = xbasis(i+1) - xbasis(i-1);
        betaY = ybasis(j+1) - ybasis(j-1);
        gradX = (pot(i+1,j)-pot(i-1,j))/betaX;
        gradY = (pot(i,j+1)-pot(i,j-1))/(betaY*xbasis(i));
        
%         velr(i,j) = gradX;
%         velp(i,j) = gradY;
        
        velr(i,j) = gradY./xbasis(i);
        velp(i,j) = -gradX;
        
        velx(i,j) = velr(i,j)*cos(ybasis(j)) - velp(i,j)*sin(ybasis(j));
        vely(i,j) = velr(i,j)*sin(ybasis(j)) + velp(i,j)*cos(ybasis(j));
        
        
    end
end

% Compute divergence

% X loop
for i=2:nX+1

    % Y loop
    for j=2:nY+1
        
        gradX = (dens(i+1,j)*velr(i+1,j)-dens(i-1,j)*velr(i-1,j))/betaX;
        gradY = (dens(i,j+1)*velp(i,j+1)-dens(i,j-1)*velp(i,j-1))/(betaY*xbasis(i));
        
        divV(i,j) = gradX + gradY + 2.0*dens(i,j)*velr(i,j)/xbasis(i) + cot(ybasis(j))*dens(i,j)*velp(i,j)/xbasis(i);
        
        
    end
end
        


%% Plot
figure
hold on;
xlim([min(X(:)),max(X(:))]);
ylim([min(Y(:)),max(Y(:))]);

ylabel('y');
xlabel('x');
modifyaxis_paper

% [gradX,gradY] = gradient(pot',hVector(1),hVector(2));
% [X,Y] = meshgrid(xbasis,ybasis);

% velx = gradY./R;
% vely = -gradX;

% velx = gradX;
% vely = gradY./R;
% 
V = sqrt(velx'.^2 + vely'.^2);
% 
% 
% [divX, ~] = gradient(dens.*velx',hVector(1));
% [~, divY] = gradient(dens.*vely',hVector(2));
% 
% div = divX' + 2*velx./R + divY'./R + cot(PHI)./R.*vely;
% 
% % div = divX + divY;
% min(div(:))
% max(div(:))

% [velX,velY] = curl(X,Y,
% 
% contour(X,Y,pot',30)

h = pcolor(X,Y,pot');
set(h,'edgealpha',0)
set(h,'facealpha',1)

la = V>=0.0*max(V(:));
la(1:3:end) = false;
VX = velx';
VY = vely';
quiver(X(la),Y(la),VX(la),VY(la),2)
% quiver(X',Y',velx,vely,1)
axis equal
            
            

