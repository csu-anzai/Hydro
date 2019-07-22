clear all; close all; clc;

nX = 200;
nY = 100;

pot = zeros(nX+2,nY+2);
source = zeros(nX+2,nY+2);


xbnds = [0,2];
ybnds = [0,1];
xbasis = linspace(xbnds(1),xbnds(2),nX+2);
ybasis = linspace(ybnds(1),ybnds(2),nY+2);

hVector = [xbasis(2)-xbasis(1);
           ybasis(2)-ybasis(1)];


[X,Y] = meshgrid(xbasis,ybasis);
Rsqr = (X-0.5).^2 + (Y-1).^2;
       
% Egg crate sources
source(:,:) = sin(pi*X').*sin(pi*Y');       
       
% Random Gaussians
source(:,:) = 0.0;
nGauss = 1000;
offset = 0.3;
for i=1:nGauss
    pos = [xbnds(1)+offset, ybnds(1)+offset] + rand(1,2).*[diff(xbnds)-2*offset,diff(ybnds)-2*offset];
    Rsqr = (X'-pos(1)).^2 + (Y'-pos(2)).^2;
    amplitude = (-1)^i;
    source = source + amplitude*exp(-1000*Rsqr);    
end

source(:,1) = 0.0;
source(:,end) = 0.0;
source(1,:) = 0.0;
source(end,:) = 0.0;
       
       
% Dirichlet of pot=1 for all x at the upper boundary
% pot(:,end) = 1.0;

iterMax = 100;

errorMax = 1e99;
errorThresh = 1e-7;

sorCoeff = 0.9;

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
            
                        
            % Compute the new Tij
            
            
            derivTerms = 0.0;
            derivTerms = derivTerms + gammaYSqr*(pot(i-1,j)+pot(i+1,j)+rX*(pot(i+1,j)-pot(i-1,j)));
            derivTerms = derivTerms + gammaXSqr*(pot(i,j-1)+pot(i,j+1)+rY*(pot(i,j+1)-pot(i,j-1)));
            
            sourceTerms = gammaXSqr*gammaYSqr*source(i,j);
            
            denom = 2.0*(gammaXSqr + gammaYSqr);
            
            
            potIJnew = (derivTerms + sourceTerms)/denom;
            
            error = abs((pot(i,j)-potIJnew));
            if(error>errorMax)
                errorMax = error;
            end
            
            pot(i,j) = (1-sorCoeff)*pot(i,j) + sorCoeff*potIJnew;
            
        end
    end
    
    fprintf('errorMax: %3.3E\n',errorMax);
    
    if(errorMax<errorThresh)
        fprintf('Converged!\n');
        break;
    end
    
end
            
            
%%
figure
hold on;
xlim(xbnds);
ylim(ybnds);
ylabel('y');
xlabel('x');
modifyaxis_paper

[gradX,gradY] = gradient(pot',hVector(1),hVector(2));
% [X,Y] = meshgrid(xbasis,ybasis);

velx = gradY;
vely = -gradX;

V = sqrt(velx.^2 + vely.^2);

% velx = gradX;
% vely = gradY;

[divX, ~] = gradient(velx,hVector(1));
[~, divY] = gradient(vely,hVector(2));
div = divX + divY;
min(div(:))
max(div(:))

% [velX,velY] = curl(X,Y,

contour(X,Y,pot',30)

la = V>0.05*max(V(:));
quiver(X(la),Y(la),velx(la),vely(la),2)

            
            

