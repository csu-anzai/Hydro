function res = integrateMonomialOverContour(coordSys, vertices, planeNormal,K,a,b,c)


% Determine the number of vertices
nVerts = size(vertices,1);

% Ensure plane normal is normalized
planeNormal = planeNormal/norm(planeNormal);


% Determine which plane to project onto based on the maximum normal
% component
% maxloc(|planeNormal|)=3 => AlphaBeta  => projectionPlane = 3
% maxloc(|planeNormal|)=2 => AlphaGamma => projectionPlane = 2
% maxloc(|planeNormal|)=1 => BetaGamma  => projectionPlane = 1
[~,ind] = max(abs(planeNormal));
if(ind==1)
    projectionPlane = 1;
elseif(ind==2)
    projectionPlane = 2;    
else
    projectionPlane = 3;
end

% Get the minimum component in the projection plane over all polygon
% vertices
% minComponent = min(vertices(:,projectionPlane));

% Project all vertices onto the projection plane at an offset of
% minComponent
vertsProj = vertices;
% vertsProj(:,projectionPlane) = minComponent;
% vertsProj

alphaInd = 1;
betaInd = 2;
gammaInd = 3;

ABplane = 3;
AGplane = 2;
BGplane = 1;

CS_CART   = 1;
CS_CYL    = 2;
CS_SPHERE = 3;

useSinIntegral = false;


currentCS = 0;
if(strcmp(coordSys,'cartesian'))
    currentCS = CS_CART;
elseif(strcmp(coordSys,'cylindrical'))
    currentCS = CS_CYL;
elseif(strcmp(coordSys,'spherical'))
    currentCS = CS_SPHERE;
else    
    fprintf(2,'[integrateContour] Coordinate system \"%s\" not supported\n',coordSys);
    return;
end

% Setup cartesian integral parameters
if(currentCS == CS_CART)
    
    % BetaGamma plane
    if(projectionPlane==BGplane)
        differentialInd = gammaInd;
        p = a;
        q = b+1;
        r = c;
        constFactor = K/(b+1);
        fprintf('Choosing BG plane\n\tDifferentialInd=%d\n\tp=%d\n\tq=%d\n\tr=%d\n\tConstant factor: %f\n\n',differentialInd,p,q,r,constFactor);
    % AlphaGamma plane
    elseif(projectionPlane==AGplane)
        differentialInd = alphaInd;
        p = a;
        q = b;
        r = c+1;
        constFactor = K/(c+1);
        fprintf('Choosing AG plane\n\tDifferentialInd=%d\n\tp=%d\n\tq=%d\n\tr=%d\n\tConstant factor: %f\n\n',differentialInd,p,q,r,constFactor);
    % AlphaBeta plane
    else
        differentialInd = betaInd;
        p = a+1;
        q = b;
        r = c;
        constFactor = K/(a+1);
        fprintf('Choosing AB plane\n\tDifferentialInd=%d\n\tp=%d\n\tq=%d\n\tr=%d\n\tConstant factor: %f\n\n',differentialInd,p,q,r,constFactor);
    end
    
    

% Setup cylindrical integral parameters
elseif(currentCS == CS_CYL)
    
    
    % BetaGamma plane
    if(projectionPlane==BGplane)
        differentialInd = gammaInd;
        p = a+1;
        q = b+1;
        r = c;
        constFactor = K/(b+1);
        
        fprintf('Choosing BG plane\n\tDifferentialInd=%d\n\tp=%d\n\tq=%d\n\tr=%d\n\tConstant factor: %f\n\n',differentialInd,p,q,r,constFactor);
    % AlphaGamma plane
    elseif(projectionPlane==AGplane)
        differentialInd = alphaInd;
        p = a;
        q = b;
        r = c+1;
        constFactor = K/(c+1);
        fprintf('Choosing AG plane\n\tDifferentialInd=%d\n\tp=%d\n\tq=%d\n\tr=%d\n\tConstant factor: %f\n\n',differentialInd,p,q,r,constFactor);
    % AlphaBeta plane
    else
        differentialInd = betaInd;
        p = a+2;
        q = b;
        r = c;
        constFactor = K/(a+2);
        fprintf('Choosing AB plane\n\tDifferentialInd=%d\n\tp=%d\n\tq=%d\n\tr=%d\n\tConstant factor: %f\n\n',differentialInd,p,q,r,constFactor);
    end
    
% Setup spherical integral parameters
elseif(currentCS == CS_SPHERE)
    
    % BetaGamma plane
    if(projectionPlane==BGplane)
        differentialInd = gammaInd;
        xiInd = gammaInd;
        useSinIntegral = true;
        p = a+2;
        q = b+1;
        r = c;
        constFactor = K/(b+1);
%         fprintf('Choosing BG plane\n\tDifferentialInd=%d\n\tp=%d\n\tq=%d\n\tr=%d\n\tConstant factor: %f\n\n',differentialInd,p,q,r,constFactor);
    % AlphaGamma plane
    elseif(projectionPlane==AGplane)
        differentialInd = gammaInd;
        p = a+2;
        q = b;
        r = c;
        constFactor = -K/(a+2);         
%         fprintf('Choosing AG plane\n\tDifferentialInd=%d\n\tp=%d\n\tq=%d\n\tr=%d\n\tConstant factor: %f\n\n',differentialInd,p,q,r,constFactor);
    % AlphaBeta plane
    else
        differentialInd = betaInd;
        xiInd = gammaInd;
        useSinIntegral = true;
        p = a+2;
        q = b;
        r = c;
%         constFactor = K*sin(minComponent)/(a+2);
        constFactor = K/(a+2);
%         fprintf('Choosing AB plane\n\tDifferentialInd=%d\n\tp=%d\n\tq=%d\n\tr=%d\n\tConstant factor: %f\n\n',differentialInd,p,q,r,constFactor);
    end
        
else
    fprintf(2,'[integrateContour] Coordinate system \"%s\" not supported\n',coordSys);
    return;
end



integral = 0;
    
% Loop over vertices
for i=1:nVerts

    % Get start and end indices for each edge
    indStart = rem(i-1,nVerts)+1;
    indEnd = rem(i,nVerts)+1;

    % Evaluate contour integral along this line segment
    % evaluateSimpleLineIntegral(muStart,muEnd,etaStart,etaEnd,p,q)
%     if(projectionPlane==BGplane && currentCS==CS_SPHERE)
    if(useSinIntegral)
        tmp = evaluateComplicatedLineIntegralSin(vertsProj(indStart,alphaInd),  vertsProj(indEnd,alphaInd),  ...
                                                 vertsProj(indStart,betaInd), vertsProj(indEnd,betaInd), ...
                                                 vertsProj(indStart,gammaInd), vertsProj(indEnd,gammaInd), ...
                                                 vertsProj(indStart,xiInd), vertsProj(indEnd,xiInd), ...
                                                 vertsProj(indStart,differentialInd), vertsProj(indEnd,differentialInd), ...
                                                 p, q, r);
    else
        tmp = evaluateSimpleLineIntegral(vertsProj(indStart,alphaInd),  vertsProj(indEnd,alphaInd),  ...
                                         vertsProj(indStart,betaInd), vertsProj(indEnd,betaInd), ...
                                         vertsProj(indStart,gammaInd), vertsProj(indEnd,gammaInd), ...
                                         vertsProj(indStart,differentialInd), vertsProj(indEnd,differentialInd), ...
                                         p, q, r);
    end
    
    integral = integral + tmp;

end

% Unproject resulting integral
% http://softsurfer.com/Archive/algorithm_0101/algorithm_0101.htm#3D Project to 2D
areaScalingFactor = 1/(planeNormal(projectionPlane)/norm(planeNormal));
% areaScalingFactor = 1;
res = constFactor*integral*areaScalingFactor;
return;




