function [gradX,gradY] = gradient2(scalar, xbasis, ybasis, nX,nY, coordSysIn)

coordSys = '';
if(nargin==5)
    coordSys = 'cartesian';
else
    coordSys = coordSysIn;
end

gradX = zeros(size(scalar));
gradY = zeros(size(scalar));

% Gradient on the interior
for i=2:nX+1
    for j=2:nY+1
        
        
        betaX = xbasis(i+1) - xbasis(i-1);
        betaY = ybasis(j+1) - ybasis(j-1);
        
        
        
        gradX(i,j) = (scalar(i+1,j) - scalar(i-1,j))/betaX;
        gradY(i,j) = (scalar(i,j+1) - scalar(i,j-1))/betaY;
        
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
        gradX(i,j) =  (-scalar(i+2,j) + 4.0*scalar(i+1,j) - 3.0*scalar(i,j))/betaX;
        gradY(i,j) =  (scalar(i,j+1) - scalar(i,j-1))/betaY;
    end
end


% Right boundary
for i=nX+2:nX+2
    for j=2:nY+1
        
        % Compute the gradients
        betaX = xbasis(i) - xbasis(i-2);
        betaY = ybasis(j+1) - ybasis(j-1);
        
        % Compute the source terms
        gradX(i,j) =  (scalar(i-2,j) - 4.0*scalar(i-1,j) + 3.0*scalar(i,j))/betaX;
        gradY(i,j) =  (scalar(i,j+1) - scalar(i,j-1))/betaY;
    end
end


% Bottom boundary
for i=2:nX+1
    for j=1:1
        
        % Compute the gradients
        betaX = xbasis(i+1) - xbasis(i-1);
        betaY = ybasis(j+2) - ybasis(j);
        
        % Compute the source terms
        gradX(i,j) =  (scalar(i+1,j) - scalar(i-1,j))/betaX;
        gradY(i,j) =  (-scalar(i,j+2) + 4.0*scalar(i,j+1) - 3.0*scalar(i,j))/betaY;
    end
end


% Top boundary
for i=2:nX+1
    for j=nY+2:nY+2
        
        % Compute the gradients
        betaX = xbasis(i+1) - xbasis(i-1);
        betaY = ybasis(j) - ybasis(j-2);
        
        % Compute the source terms
        gradX(i,j) =  (scalar(i+1,j) - scalar(i-1,j))/betaX;
        gradY(i,j) =  (scalar(i,j-2) - 4.0*scalar(i,j-1) + 3.0*scalar(i,j))/betaY;
    end
end

% Bottom-left cell
betaX = xbasis(3) - xbasis(1);
betaY = ybasis(3) - ybasis(1);

gradX(1,1) =  (-scalar(3,1) + 4.0*scalar(2,1) - 3.0*scalar(1,1))/betaX;
gradY(1,1) =  (-scalar(1,3) + 4.0*scalar(1,2) - 3.0*scalar(1,1))/betaY;


% Top-left cell
betaX = xbasis(3) - xbasis(1);
betaY = ybasis(end) - ybasis(end-2);

gradX(1,end) =  (-scalar(3,end) + 4.0*scalar(2,end) - 3.0*scalar(1,end))/betaX;
gradY(1,end) =  (scalar(1,end-2) - 4.0*scalar(1,end-1) + 3.0*scalar(1,end))/betaY;


% Bottom-right cell
betaX = xbasis(end) - xbasis(end-2);
betaY = ybasis(3) - ybasis(1);

gradX(end,1) =  (scalar(end-2,1) - 4.0*scalar(end-1,1) + 3.0*scalar(end,1))/betaX;
gradY(end,1) =  (-scalar(end,3) + 4.0*scalar(end,2) - 3.0*scalar(end,1))/betaY;


% Top-right cell
betaX = xbasis(end) - xbasis(end-2);
betaY = ybasis(end) - ybasis(end-2);

gradX(end,end) =  (scalar(end-2,end) - 4.0*scalar(end-1,end) + 3.0*scalar(end,end))/betaX;
gradY(end,end) =  (scalar(end,end-2) - 4.0*scalar(end,end-1) + 3.0*scalar(end,end))/betaY;

% Do coordinate system corrections
switch(coordSys)
    case('cartesian')
        % Do nothing
    case('spherical')
        % Divide the gradient in Y by the X coordinate
        for i=2:nX+1
            for j=2:nY+1
                gradY(i,j) = gradY(i,j)/xbasis(i);
            end
        end
                
    otherwise
        fprintf(2,'[gradient2] ERROR: Coordinate system \"%s\" not supported!\n',coordSys);
end
        







