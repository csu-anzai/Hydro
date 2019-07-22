function div = divergence2(vX, vY, xbasis, ybasis, nX,nY, coordSysIn)

coordSys = '';
if(nargin==6)
    coordSys = 'cartesian';
else
    coordSys = coordSysIn;
end

div = zeros(size(vX));

% Divergence on the interior
for i=2:nX+1
    for j=2:nY+1
        
        
        betaX = xbasis(i+1) - xbasis(i-1);
        betaY = ybasis(j+1) - ybasis(j-1);
        
        
        
        
        switch(coordSys)
            case('cartesian')
                gradX = (vX(i+1,j) - vX(i-1,j))/betaX;
                gradY = (vY(i,j+1) - vY(i,j-1))/betaY;
                div(i,j) = gradX + gradY;
            case('spherical')
                r    = xbasis(i);
                rSqr = r;
                rP1  = xbasis(i+1);
                rM1  = xbasis(i-1);
                p    = ybasis(j);
                pP1  = ybasis(j+1);
                pM1  = ybasis(j-1);
                
                gradX = (1/rSqr)*(rP1*rP1*vX(i+1,j) - rM1*rM1*vX(i-1,j))/betaX;
                gradY = (1/(r*sin(p)))*(sin(pP1)*vY(i,j+1) - sin(pM1)*vY(i,j-1))/betaY;
                
                div(i,j) = gradX + gradY;
                
            otherwise
                fprintf(2,'[gradient2] ERROR: Coordinate system \"%s\" not supported!\n',coordSys);
        end
        
    end
end


% Compute the gradient on the boundary via forward/backward differences

% Left boundary
for i=1:1
    for j=2:nY+1
        
        % Compute the gradients
        betaX = xbasis(i+2) - xbasis(i);
        betaY = ybasis(j+1) - ybasis(j-1);
        
        % Compute the source terms
        gradX =  (-vX(i+2,j) + 4.0*vX(i+1,j) - 3.0*vX(i,j))/betaX;
        gradY =  (vY(i,j+1) - vY(i,j-1))/betaY;
        
        switch(coordSys)
            case('cartesian')
                div(i,j) = gradX + gradY;
            case('spherical')
                div(i,j) = gradX + (gradY + 2.0*vX(i,j) + cot(ybasis(j))*vY(i,j))/xbasis(i);
            otherwise
                fprintf(2,'[gradient2] ERROR: Coordinate system \"%s\" not supported!\n',coordSys);
        end
    end
end


% Right boundary
for i=nX+2:nX+2
    for j=2:nY+1
        
        % Compute the gradients
        betaX = xbasis(i) - xbasis(i-2);
        betaY = ybasis(j+1) - ybasis(j-1);
        
        % Compute the source terms
        gradX =  (vX(i-2,j) - 4.0*vX(i-1,j) + 3.0*vX(i,j))/betaX;
        gradY =  (vY(i,j+1) - vY(i,j-1))/betaY;
        switch(coordSys)
            case('cartesian')
                div(i,j) = gradX + gradY;
            case('spherical')
                div(i,j) = gradX + (gradY + 2.0*vX(i,j) + cot(ybasis(j))*vY(i,j))/xbasis(i);
            otherwise
                fprintf(2,'[gradient2] ERROR: Coordinate system \"%s\" not supported!\n',coordSys);
        end
    end
end


% Bottom boundary
for i=2:nX+1
    for j=1:1
        
        % Compute the gradients
        betaX = xbasis(i+1) - xbasis(i-1);
        betaY = ybasis(j+2) - ybasis(j);
        
        % Compute the source terms
        gradX =  (vX(i+1,j) - vX(i-1,j))/betaX;
        gradY =  (-vY(i,j+2) + 4.0*vY(i,j+1) - 3.0*vY(i,j))/betaY;
        switch(coordSys)
            case('cartesian')
                div(i,j) = gradX + gradY;
            case('spherical')
                div(i,j) = gradX + (gradY + 2.0*vX(i,j) + cot(ybasis(j))*vY(i,j))/xbasis(i);
            otherwise
                fprintf(2,'[gradient2] ERROR: Coordinate system \"%s\" not supported!\n',coordSys);
        end
    end
end


% Top boundary
for i=2:nX+1
    for j=nY+2:nY+2
        
        % Compute the gradients
        betaX = xbasis(i+1) - xbasis(i-1);
        betaY = ybasis(j) - ybasis(j-2);
        
        % Compute the source terms
        gradX =  (vX(i+1,j) - vX(i-1,j))/betaX;
        gradY =  (vY(i,j-2) - 4.0*vY(i,j-1) + 3.0*vY(i,j))/betaY;
        switch(coordSys)
            case('cartesian')
                div(i,j) = gradX + gradY;
            case('spherical')
                div(i,j) = gradX + (gradY + 2.0*vX(i,j) + cot(ybasis(j))*vY(i,j))/xbasis(i);
            otherwise
                fprintf(2,'[gradient2] ERROR: Coordinate system \"%s\" not supported!\n',coordSys);
        end
    end
end

% Bottom-left cell
betaX = xbasis(3) - xbasis(1);
betaY = ybasis(3) - ybasis(1);

gradX =  (-vX(3,1) + 4.0*vX(2,1) - 3.0*vX(1,1))/betaX;
gradY =  (-vY(1,3) + 4.0*vY(1,2) - 3.0*vY(1,1))/betaY;
switch(coordSys)
    case('cartesian')
        div(1,1) = gradX + gradY;
    case('spherical')
        div(1,1) = gradX + (gradY + 2.0*vX(1,1) + cot(ybasis(1))*vY(1,1))/xbasis(1);
    otherwise
        fprintf(2,'[gradient2] ERROR: Coordinate system \"%s\" not supported!\n',coordSys);
end


% Top-left cell
betaX = xbasis(3) - xbasis(1);
betaY = ybasis(end) - ybasis(end-2);

gradX =  (-vX(3,end) + 4.0*vX(2,end) - 3.0*vX(1,end))/betaX;
gradY =  (vY(1,end-2) - 4.0*vY(1,end-1) + 3.0*vY(1,end))/betaY;
switch(coordSys)
    case('cartesian')
        div(1,end) = gradX + gradY;
    case('spherical')
        div(1,end) = gradX + (gradY + 2.0*vX(1,end) + cot(ybasis(end))*vY(1,end))/xbasis(1);
    otherwise
        fprintf(2,'[gradient2] ERROR: Coordinate system \"%s\" not supported!\n',coordSys);
end

% Bottom-right cell
betaX = xbasis(end) - xbasis(end-2);
betaY = ybasis(3) - ybasis(1);

gradX =  (vX(end-2,1) - 4.0*vX(end-1,1) + 3.0*vX(end,1))/betaX;
gradY =  (-vY(end,3) + 4.0*vY(end,2) - 3.0*vY(end,1))/betaY;
switch(coordSys)
    case('cartesian')
        div(end,1) = gradX + gradY;
    case('spherical')
        div(end,1) = gradX + (gradY + 2.0*vX(end,1) + cot(ybasis(1))*vY(end,1))/xbasis(end);
    otherwise
        fprintf(2,'[gradient2] ERROR: Coordinate system \"%s\" not supported!\n',coordSys);
end

% Top-right cell
betaX = xbasis(end) - xbasis(end-2);
betaY = ybasis(end) - ybasis(end-2);

gradX =  (vX(end-2,end) - 4.0*vX(end-1,end) + 3.0*vX(end,end))/betaX;
gradY =  (vY(end,end-2) - 4.0*vY(end,end-1) + 3.0*vY(end,end))/betaY;
switch(coordSys)
    case('cartesian')
        div(end,end) = gradX + gradY;
    case('spherical')
        div(end,end) = gradX + (gradY + 2.0*vX(end,end) + cot(ybasis(end))*vY(end,end))/xbasis(end);
    otherwise
        fprintf(2,'[gradient2] ERROR: Coordinate system \"%s\" not supported!\n',coordSys);
end

        







