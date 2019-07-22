function potential = solvePoisson(potentialInitial, RHS,xbasis,ybasis,nX,nY,bcType,errorThresh,iterMax)

potential = potentialInitial;

% errorThresh = 1e-5;
sorCoeff = 0.9;
% iterMax = 10000;

errorCoords = [0,0];

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
            
            dxSqr = (xbasis(i)-xbasis(i-1))*(xbasis(i+1)-xbasis(i));
            dySqr = (ybasis(i)-ybasis(i-1))*(ybasis(i+1)-ybasis(i));
            
            potLapla = 0.0;
            potLapla = potLapla + dySqr*(potential(i+1,j)+potential(i-1,j));
            potLapla = potLapla + dxSqr*(potential(i,j+1)+potential(i,j-1));
            
            
            % source terms
            source = dxSqr*dySqr*RHS(i,j);
            
            % prefactor on the IJ term
            denom = -2*(dxSqr+dySqr);

            
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
    
    % Y loop
    for j=nY+1:-1:2
        % X loop
        for i=nX+1:-1:2
        
            % http://ocw.usu.edu/civil_and_environmental_engineering/numerical_methods_in_civil_engineering/LaplaceDirichletTemperature.pdf
            % Compute hx^2, hy^2, and beta^2, and the denominator for
            % the old terms
            
            dxSqr = (xbasis(i)-xbasis(i-1))*(xbasis(i+1)-xbasis(i));
            dySqr = (ybasis(i)-ybasis(i-1))*(ybasis(i+1)-ybasis(i));
            
            potLapla = 0.0;
            potLapla = potLapla + dySqr*(potential(i+1,j)+potential(i-1,j));
            potLapla = potLapla + dxSqr*(potential(i,j+1)+potential(i,j-1));
            
            
            % source terms
            source = dxSqr*dySqr*RHS(i,j);
            
            % prefactor on the IJ term
            denom = -2*(dxSqr+dySqr);

            
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
    
    switch(bcType)
        case{'dir','dirichlet'}
            % Do nothing, assume the proper BC are set in the initial data
        case{'per','periodic'}
            potential(:,end) = potential(:,2);
            potential(:,1)   = potential(:,end-1);
            potential(1,:)   = potential(end-1,:);
            potential(end,:) = potential(2,:);
        case{'neu','neumann','outflow'}
            potential(:,end) = potential(:,end-1);
            potential(:,1)   = potential(:,2); 
            potential(1,:)   = potential(2,:); 
            potential(end,:) = potential(end-1,:);
%             return;
        otherwise
            fprintf(2,'[solvePoisson] ERROR: bcType = \"%s\" not supported!\n',bcType);
            return;
    end
        
%     potential(:,end) = potential(:,end-1); % Neumann on top 
%     potential(:,1) = potential(:,2); % Neumann on bottom 
%     potential(1,:) = 0;
%     potential(end,:) = 0;
%     
    if(rem(iter,100)==0 || iter == 1)
        fprintf('iteration %d\n',iter);
        fprintf('errorMax: %3.3E (%d , %d)\n',errorMax,errorCoords);
    end
    
    if(errorMax<errorThresh)
        fprintf('Converged! %d %E (%E)\n',iter, errorMax, errorThresh);
        break;
    end
    
end