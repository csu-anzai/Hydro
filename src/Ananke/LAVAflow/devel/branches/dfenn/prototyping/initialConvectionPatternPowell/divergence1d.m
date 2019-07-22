function divV = divergence1d(V,xbasis,nX,coordSysIn)

divV = zeros(size(V));



coordSys = '';
if(nargin==4)
    coordSys = 'cartesian';
else
    coordSys = coordSysIn;
end

div = zeros(size(V));

% Divergence on the interior
for i=2:nX-1
        
        
    betaX = xbasis(i+1) - xbasis(i-1);      



    switch(coordSys)
        case('cartesian')
            gradX = (V(i+1) - V(i-1))/betaX;
            divV(i) = gradX;
        case('spherical')
            r    = xbasis(i);
            rSqr = r;
            rP1  = xbasis(i+1);
            rM1  = xbasis(i-1);

            gradX = (1/rSqr)*(rP1*rP1*V(i+1) - rM1*rM1*V(i-1))/betaX;

            divV(i) = gradX;

        otherwise
            fprintf(2,'[gradient2] ERROR: Coordinate system \"%s\" not supported!\n',coordSys);
    end

end



% % Divergence on the interior
% for i=2:nX-1
%     betaX = xbasis(i+1) - xbasis(i-1);
%     divV(i) = (V(i+1) - V(i-1))/betaX;
% end

% Divergence on the left
switch(coordSys)
    case('cartesian')
        betaX = xbasis(3) - xbasis(1);
        divV(1) =  (-V(3) + 4.0*V(2) - 3.0*V(1))/betaX;
    case('spherical')
        r = xbasis(1);
        rP1 = xbasis(2);
        rSqr = r*r;
        betaX = xbasis(2)-xbasis(1);
        divV(1) = (rP1*rP1*V(2) - r*r*V(1))/(rSqr*betaX);
end


% Divergence on the right
switch(coordSys)
    case('cartesian')
        betaX = xbasis(end) - xbasis(end-2);
        divV(end) =  (V(end-2) - 4.0*V(end-1) + 3.0*V(end))/betaX;
    case('spherical')
        r = xbasis(end);
        rM1 = xbasis(end-1);
        rSqr = r*r;
        betaX = xbasis(end)-xbasis(end-1);
        divV(end) = (r*r*V(end) - rM1*rM1*V(end-1))/(rSqr*betaX);
end