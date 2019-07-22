function [potential,errorMax] = solvePoissonLargeStencil(potentialInitial, RHS,xbasis,ybasis,nX,nY,bcType,errorThresh,iterMax, coordSysIn)

coordSys = '';
if(nargin==9)
    coordSys = 'cartesian';
else
    coordSys = coordSysIn;
end

potential = potentialInitial;

% errorThresh = 1e-5;
sorCoeff = 0.9;
% iterMax = 10000;

errorCoords = [0,0];

% Iteration Loop
for iter=1:iterMax
    
%     fprintf('iteration %d\n',iter);
    errorMax = -1;
    
        switch(bcType)
        case{'dir','dirichlet'}
            % Do nothing, assume the proper BC are set in the initial data
        case{'per','periodic'}
            potential(:,end) = potential(:,4);
            potential(:,end-1) = potential(:,3);
            
            potential(:,1)   = potential(:,end-3);
            potential(:,2)   = potential(:,end-2);
            
            potential(1,:)   = potential(end-3,:);
            potential(2,:)   = potential(end-2,:);
            
            potential(end,:) = potential(4,:);
            potential(end-1,:) = potential(3,:);
                        
        case{'neu','neumann','outflow'}
            potential(:,end) = potential(:,end-2);
            potential(:,end-1) = potential(:,end-2);
            potential(:,1)   = potential(:,3); 
            potential(:,2)   = potential(:,3); 
            potential(1,:)   = potential(3,:); 
            potential(2,:)   = potential(3,:); 
            potential(end,:) = potential(end-2,:);
            potential(end-1,:) = potential(end-2,:);
%             return;
        case{'solenoidal'}
            % This case assumes that the normal dotted with curl(vector
            % potential) is zero
            %
            % This reduces the Hodge decomposition to a Neumann BC
            % dot(n,v) = dot(n, grad(potential)) = d(potential)/d(n) on the boundary
            %
            % We assume that the dot(n,v) term is stored in the boundary
            % cells of the initial potential field
            
            
            % Bottom boundary (r=rmin)
            potPrime2 = potentialInitial(1,:);
            potPrime1 = potentialInitial(2,:);
            a = 0.5*(potPrime1-potPrime2)/(xbasis(2)-xbasis(1));
            b = potPrime2 - a*xbasis(1);
            c = potential(3,:) - a*xbasis(3)^2 - b*xbasis(3);
            potential(1,:) = a*xbasis(1) + b*xbasis(1) + c;
            potential(2,:) = a*xbasis(2) + b*xbasis(2) + c;
            
            
            % Top boundary (r=rmax)
            potPrime2 = potentialInitial(end-1,:);
            potPrime1 = potentialInitial(end,:);
            a = 0.5*(potPrime1-potPrime2)/(xbasis(end-1)-xbasis(end));
            b = potPrime2 - a*xbasis(end);
            c = potential(end-2,:) - a*xbasis(end-2)^2 - b*xbasis(end-2);
            potential(end,:) = a*xbasis(end) + b*xbasis(end) + c;
            potential(end-1,:) = a*xbasis(end-1) + b*xbasis(end-1) + c;
            
            
            
            % Left boundary (phi=phimin)
            potPrime2 = potentialInitial(:,1);
            potPrime1 = potentialInitial(:,2);
            a = 0.5*(potPrime1-potPrime2)/(ybasis(2)-ybasis(1));
            b = potPrime2 - a*ybasis(1);
            c = potential(:,3) - a*ybasis(3)^2 - b*ybasis(3);
            potential(:,1) = a*ybasis(1) + b*ybasis(1) + c;
            potential(:,2) = a*ybasis(2) + b*ybasis(2) + c;
            
            % Right boundary (phi=phimax)
            potPrime2 = potentialInitial(:,end);
            potPrime1 = potentialInitial(:,end-1);
            a = 0.5*(potPrime1-potPrime2)/(ybasis(end-1)-ybasis(end));
            b = potPrime2 - a*ybasis(end);
            c = potential(:,end-2) - a*ybasis(end-2)^2 - b*ybasis(end-2);
            potential(:,end) = a*ybasis(end) + b*ybasis(end) + c;
            potential(:,end-1) = a*ybasis(end-1) + b*ybasis(end-1) + c;
            
        case{'supernova','toth'}
            % This case takes advice from Toth (2000) for the boundary
            % conditions for elliptic cleaning
            
            intX = 3:nX+2;
            intY = 3:nY+2;
            
%             % Set the lower boundary to fixed flow
%             potential(2,intY) = 0;
%             potential(1,intY) = potential(3,intY);

            % Set the lower boundary to open
            potential(2,intY) = 0;
            potential(1,intY) = 0;

%             % Set the lower boundary to symmetric
%             potential(1,intY) = potential(4,intY);
%             potential(2,intY) = potential(3,intY);
            
            % Set the upper boundary to open
            potential(end,intY) = 0;
            potential(end-1,intY) = 0;
            
%             % Set the upper boundary to fixed flow
%             potential(end,intY) = 0;
%             potential(end-1,intY) = potential(3,intY);
            
            % Set the left boundary to symmetric
            potential(intX,end)   = potential(intX,end-3);
            potential(intX,end-1) = potential(intX,end-2);
            
            % Set the right boundary to symmetric
            potential(intX,1) = potential(intX,4);
            potential(intX,2) = potential(intX,3);
            
            
            
        otherwise
            fprintf(2,'[solvePoisson] ERROR: bcType = \"%s\" not supported!\n',bcType);
            return;
    end
    
    % X loop
    for i=3:nX+2
        
        % Y loop
        for j=3:nY+2
            % http://ocw.usu.edu/civil_and_environmental_engineering/numerical_methods_in_civil_engineering/LaplaceDirichletTemperature.pdf
            % Compute hx^2, hy^2, and beta^2, and the denominator for
            % the old terms
            dx = xbasis(i+1)-xbasis(i);
            dy = ybasis(j+1)-ybasis(j);
            dxSqr = (xbasis(i)-xbasis(i-1))*(xbasis(i+1)-xbasis(i));
            dySqr = (ybasis(i)-ybasis(i-1))*(ybasis(i+1)-ybasis(i));
            
            potLapla = 0.0;
            source = 0.0;
            denom = 0.0;
            switch(coordSys)
                case('cartesian')
                    potLapla = potLapla + dySqr*(potential(i+2,j)+potential(i-2,j));
                    potLapla = potLapla + dxSqr*(potential(i,j+2)+potential(i,j-2));
                    
                    % source terms
                    source = source + 4.0*dxSqr*dySqr*RHS(i,j);
                    
                    % prefactor on the IJ term
                    denom = -2*(dxSqr+dySqr);
            
                case('spherical')
                    r = xbasis(i);
                    rSqr = r*r;
                    rM1 = xbasis(i-1);
                    rP1 = xbasis(i+1);
                                        
                    phi = ybasis(j);
                    phiM1 = ybasis(j-1);
                    phiP1 = ybasis(j+1);
                    
                    betaX = xbasis(i+1)-xbasis(i-1);
                    betaY = ybasis(j+1)-ybasis(j-1);
                    
%                     potLapla = potLapla + dySqr*(potential(i+2,j)+potential(i-2,j));
%                     potLapla = potLapla + dxSqr*(potential(i,j+2)+potential(i,j-2))/rSqr;
%                     potLapla = potLapla + 4.0*dxSqr*dySqr*(2.0/r)*(potential(i+1,j)-potential(i-1,j))/betaX;
%                     potLapla = potLapla + 4.0*dxSqr*dySqr*(cot(phi)/rSqr)*(potential(i,j+1)-potential(i,j-1))/betaY;
                    
%                     potLapla = potLapla + 4.0*(2.0/r)*dxSqr*dySqr*(-potential(i+2,j)+8*potential(i+1,j)-8*potential(i-1,j)+potential(i-2,j))/(12*dx);
%                     potLapla = potLapla + 4.0*dxSqr*dySqr*(cot(phi)/rSqr)*(-potential(i,j+2)+8*potential(i,j+1)-8*potential(i,j-1)+potential(i,j-2))/(12*dy);
                    

                    potLapla = potLapla + dySqr*((rP1*rP1 - rM1*rM1)/rSqr)*(potential(i+1,j)-potential(i-1,j));
                    potLapla = potLapla + dySqr*(potential(i+2,j)+potential(i-2,j));
                    potLapla = potLapla + dxSqr*((sin(phiP1)/rP1 - sin(phiM1)/rM1)/(r*sin(phi)))*(potential(i,j+1)-potential(i,j-1));
                    potLapla = potLapla + (dxSqr/rSqr)*(potential(i,j+2)+potential(i,j-2));

                    
                    % source terms
                    source = source + 4.0*dxSqr*dySqr*RHS(i,j);
                    
                    % prefactor on the IJ term
                    denom = -2*(dxSqr/rSqr + dySqr);
                    
                otherwise
                    fprintf(2,'[solvePoissonLargeStencil] ERROR: Coordinate system \"%s\" not supported!\n',coordSys);
            end
            
            

            
            if(~isreal(potLapla))
                fprintf(2,'Laplacian of potential is complex!\n');
            end
            if(~isreal(source))
                fprintf(2,'RHS is complex!\n');
            end
            
            
            potNew = (source - potLapla)/denom;
                       
            
            
            % Get the old values
            potOld = potential(i,j);
            potNew = (1-sorCoeff)*potOld + sorCoeff*potNew;
            
            
            if(~isreal(source))
                fprintf(2,'New potential is complex!\n');
            end
            
            err = abs(potNew-potOld)/abs(potNew);
            
            if(err>errorMax && ~isnan(err))
                errorMax = err;
                errorCoords = [i,j];
            end
            
            potential(i,j) = potNew;
            
        end
    end
    
%     % Y loop
%     for j=nY+2:-1:3
%         % X loop
%         for i=nX+2:-1:3
%         
%             % http://ocw.usu.edu/civil_and_environmental_engineering/numerical_methods_in_civil_engineering/LaplaceDirichletTemperature.pdf
%             % Compute hx^2, hy^2, and beta^2, and the denominator for
%             % the old terms
%             
%             dx = xbasis(i+1)-xbasis(i);
%             dy = ybasis(j+1)-ybasis(j);
%             dxSqr = (xbasis(i)-xbasis(i-1))*(xbasis(i+1)-xbasis(i));
%             dySqr = (ybasis(i)-ybasis(i-1))*(ybasis(i+1)-ybasis(i));
%             
%             potLapla = 0.0;
% %             potLapla = potLapla + dySqr*(potential(i+2,j)+potential(i-2,j));
% %             potLapla = potLapla + dxSqr*(potential(i,j+2)+potential(i,j-2));
%             switch(coordSys)
%                 case('cartesian')
%                     potLapla = potLapla + dySqr*(potential(i+2,j)+potential(i-2,j));
%                     potLapla = potLapla + dxSqr*(potential(i,j+2)+potential(i,j-2));
%                     
%                     % source terms
%                     source = source + 4.0*dxSqr*dySqr*RHS(i,j);
%                     
%                     % prefactor on the IJ term
%                     denom = -2*(dxSqr+dySqr);
%                     
%                 case('spherical')
%                     r = xbasis(i);
%                     rSqr = r*r;
%                     rM1 = xbasis(i-1);
%                     rP1 = xbasis(i+1);
%                                         
%                     phi = ybasis(j);
%                     phiM1 = ybasis(j-1);
%                     phiP1 = ybasis(j+1);
%                     
%                     betaX = xbasis(i+1)-xbasis(i-1);
%                     betaY = ybasis(j+1)-ybasis(j-1);
%                     
% %                     potLapla = potLapla + dySqr*(potential(i+2,j)+potential(i-2,j));
% %                     potLapla = potLapla + dxSqr*(potential(i,j+2)+potential(i,j-2))/rSqr;
% %                     potLapla = potLapla + 4.0*(2.0/r)*dxSqr*dySqr*(potential(i+1,j)-potential(i-1,j))/betaX;
% %                     potLapla = potLapla + 4.0*dxSqr*dySqr*(cot(phi)/rSqr)*(potential(i,j+1)-potential(i,j-1))/betaY;
% %                     
% %                     potLapla = potLapla + 4.0*(2.0/r)*dxSqr*dySqr*(-potential(i+2,j)+8*potential(i+1,j)-8*potential(i-1,j)+potential(i-2,j))/(12*dx);
% %                     potLapla = potLapla + 4.0*dxSqr*dySqr*(cot(phi)/rSqr)*(-potential(i,j+2)+8*potential(i,j+1)-8*potential(i,j-1)+potential(i,j-2))/(12*dy);
% 
% 
%                     potLapla = potLapla + dySqr*((rP1*rP1 - rM1*rM1)/rSqr)*(potential(i+1,j)-potential(i-1,j));
%                     potLapla = potLapla + dySqr*(potential(i+2,j)+potential(i-2,j));
%                     potLapla = potLapla + dxSqr*((sin(phiP1)/rP1 - sin(phiM1)/rM1)/(r*sin(phi)))*(potential(i,j+1)-potential(i,j-1));
%                     potLapla = potLapla + (dxSqr/rSqr)*(potential(i,j+2)+potential(i,j-2));
% 
%                     
%                     % source terms
%                     source = source + 4.0*dxSqr*dySqr*RHS(i,j);
%                     
%                     % prefactor on the IJ term
%                     denom = -2*(dxSqr/rSqr + dySqr);
%                     
%                 otherwise
%                     fprintf(2,'[solvePoissonLargeStencil] ERROR: Coordinate system \"%s\" not supported!\n',coordSys);
%             end
%             
% %             % source terms
% %             source = 4.0*dxSqr*dySqr*RHS(i,j);
% %             
% %             % prefactor on the IJ term
% %             denom = -2*(dxSqr+dySqr);
% 
%             
%             if(~isreal(potLapla))
%                 fprintf(2,'Laplacian of potential is complex!\n');
%             end
%             if(~isreal(source))
%                 fprintf(2,'RHS is complex!\n');
%             end
%             
%             
%             potNew = (source - potLapla)/denom;
%                        
%             
%             
%             % Get the old values
%             potOld = potential(i,j);
%             potNew = (1-sorCoeff)*potOld + sorCoeff*potNew;
%             
%             
%             if(~isreal(source))
%                 fprintf(2,'New potential is complex!\n');
%             end
%             
%             err = abs(potNew-potOld)/abs(potNew);
%             
%             if(err>errorMax && ~isnan(err))
%                 errorMax = err;
%                 errorCoords = [i,j];
%             end
%             
%             potential(i,j) = potNew;
%             
%         end
%     end
    
    % Set boundary conditions
%     potential(:,end) = potential(:,end-1); % Neumann on top 
%     potential(:,1) = potential(:,2); % Neumann on bottom 
%     divPot(1,:) = divPot(2,:); % Neumann on left 
%     divPot(end,:) = divPot(end-1,:); % Neumann on right 
%     divPot(:,end) = 0;
%     divPot(:,1) = 0;
%     potential(1,:) = 0;
%     potential(end,:) = 0;
    
    
    
%     divPot(:,end) = divPot(:,end-1); % Neumann on top 
%     divPot(:,1) = divPot(:,2); % Neumann on bottom 
%     divPot(1,:) = divPot(2,:) + (divPot(3,:)-divPot(2,:))/(xbasis(3)-xbasis(2))*(xbasis(1)-xbasis(2));
%     divPot(end,:) = divPot(end-2,:) + (divPot(end-1,:)-divPot(end-2,:))/(xbasis(end-1)-xbasis(end-1))*(xbasis(end)-xbasis(end-2));
    

        
%     potential(:,end) = potential(:,end-1); % Neumann on top 
%     potential(:,1) = potential(:,2); % Neumann on bottom 
%     potential(1,:) = 0;
%     potential(end,:) = 0;
%     
    if(rem(iter,100)==0)% || iter == 1)
        fprintf('iteration %d\n',iter);
        fprintf('errorMax: %3.3E (%d , %d)\n',errorMax,errorCoords);
    end
    
    if(errorMax<errorThresh)
        fprintf('Converged! %d %E (%E)\n',iter, errorMax, errorThresh);
        break;
    end
    
end