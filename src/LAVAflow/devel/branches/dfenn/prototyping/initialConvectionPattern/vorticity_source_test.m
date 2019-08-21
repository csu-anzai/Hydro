clear all; close all; clc;

nX = 100;
nY = 100;

velx = ones(nX+2,nY+2);
vely = ones(nX+2,nY+2);
vorticity = zeros(nX+2,nY+2);



% vorticity(20:40,20:40) = 1.0;





xbnds = [0,1];
ybnds = [0,2];
xbasis = linspace(xbnds(1),xbnds(2),nX+2);
ybasis = linspace(ybnds(1),ybnds(2),nY+2);

hVector = [xbasis(2)-xbasis(1);
           ybasis(2)-ybasis(1)];

[X,Y] = meshgrid(xbasis,ybasis);
Rsqr = (X-0.5).^2 + (Y-1).^2;
       
% vorticity(:,:) = 1000*sin(2*pi*X).*sin(2*pi*Y).*exp(-100*Rsqr);
vorticity(50:60,50:60) = 1000;
       
       
       
% Dirichlet of pot=1 for all x at the upper boundary
% pot(:,end) = 1.0;

iterMax = 400;

errorMax = 1e99;
errorThresh = 1e-4;

sorCoeff = 0.05;

isForward = true;
% Iteration Loop
for iter=1:iterMax
    
    fprintf('iteration %d\n',iter);
    errorMaxX = -1;
    errorMaxY = -1;
    
    isForward = rem(iter,2)==0;
    
    offset = 1;
    if(~isForward)
        offset = -1;
    end
    
    % X loop
    for i=2:nX+1
        
        % Y loop
        for j=2:nY+1
            
            % Compute spatial terms
%             dxSqr = (xbasis(i+1)-xbasis(i-1))^2;
%             dySqr = (ybasis(j+1)-ybasis(j-1))^2;
            
            RHS = zeros(2,1);
            
            if(~isForward)
                dx = xbasis(i)-xbasis(i-1);
                dy = ybasis(j)-ybasis(j-1);

                % Set up matrix to solve for velx and vely
                % This comes from the equation for vorticity and 
                % the incompressible constraint
                A = [-dx, dy; dy, dx];

                % Set up RHS
                RHS(1) = dx*dy*vorticity(i,j) - dx*velx(i,j-1) + dy*vely(i-1,j);
                RHS(2) = dy*velx(i-1,j) + dx*vely(i,j-1);

            else
                dx = xbasis(i+1)-xbasis(i);
                dy = ybasis(j+1)-ybasis(j);

                % Set up matrix to solve for velx and vely
                % This comes from the equation for vorticity and 
                % the incompressible constraint
                A = [dx, -dy; dy, dx];

                % Set up RHS
                RHS(1) = dx*dy*vorticity(i,j) + dx*velx(i,j+1) - dy*vely(i+1,j);
                RHS(2) = dy*velx(i+1,j) + dx*vely(i,j+1);
            end
            
            
            % Solve for velxIJ and velyIJ
            res = A\RHS;
%             res = linsolve(A,RHS);

            velxIJ = res(1);
            velyIJ = res(2);


%             vely(i,j) = velyIJ;
            
            errorX = abs(velx(i,j)-velxIJ);%/(abs(velx(i,j)) + abs(velxIJ));
            if(errorX>errorMaxX)
                errorMaxX = errorX;
            end
            errorY = abs(vely(i,j)-velyIJ);%/(abs(vely(i,j)) + abs(velyIJ));
            if(errorY>errorMaxY)
                errorMaxY = errorY;
            end

            
            velx(i,j) = (1-sorCoeff)*velx(i,j) + sorCoeff*velxIJ;
            vely(i,j) = (1-sorCoeff)*vely(i,j) + sorCoeff*velyIJ;
            
            % Flip isForward
%             isForward = ~isForward;
            
%             potIJnew = (derivTerms + sourceTerms)/denom;
%             
%             error = abs((pot(i,j)-potIJnew));
%             if(error>errorMax)
%                 errorMax = error;
%             end
%             
%             pot(i,j) = potIJnew;
            
        end
    end
    
    fprintf('errorMaxX: %3.3E\n',errorMaxX);
    fprintf('errorMaxY: %3.3E\n',errorMaxY);
    
    if(errorMaxX<errorThresh && errorMaxY<errorThresh)
        fprintf('Converged!\n');
        break;
    end
    
end
            
            
%%
figure
hold on;

% [gradX,gradY] = gradient(pot,hVector(1),hVector(2));
[X,Y] = meshgrid(xbasis,ybasis);

% [velX,velY] = curl(X,Y,

contour(X,Y,vorticity)
quiver(X,Y,velx,vely,1)

%% Compute the divergence of the velocity field
[divX, ~] = gradient(velx,hVector(1));
[~, divY] = gradient(vely,hVector(2));

div = divX + divY;

min(div(:))
max(div(:))
            

            
            

