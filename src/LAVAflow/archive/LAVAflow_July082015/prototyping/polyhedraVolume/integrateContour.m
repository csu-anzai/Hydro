function res = integrateContour(coordSys, vertices, planeNormal)


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
minComponent = min(vertices(:,projectionPlane));

% Project all vertices onto the projection plane at an offset of
% minComponent
vertsProj = vertices;
vertsProj(:,projectionPlane) = minComponent;


% Do cartesian stuff
if(strcmp(coordSys,'cartesian'))
    
    % BetaGamma plane
    if(projectionPlane==1)
        muInd = 2;
        etaInd = 3;
        differentialInd = 3;
        p=1;
        q = 0;
        constFactor = 1;
        fprintf('Choosing BG plane\n\tmuInd=%d\n\tetaInd=%d\n\tp=%d\n\tq=%d\n\tConstant factor: %f\n\n',muInd,etaInd,p,q,constFactor);
    % AlphaGamma plane
    elseif(projectionPlane==2)
        muInd = 3;
        etaInd = 1;
        differentialInd = 1;
        p=1;
        q=0;
        constFactor = 1;
        fprintf('Choosing AG plane\n\tmuInd=%d\n\tetaInd=%d\n\tp=%d\n\tq=%d\n\tConstant factor: %f\n\n',muInd,etaInd,p,q,constFactor);
    % AlphaBeta plane
    else
        muInd = 1;
        etaInd = 2;
        differentialInd = 2;
        p=1;
        q=0;
        constFactor = 1;
        fprintf('Choosing AB plane\n\tmuInd=%d\n\tetaInd=%d\n\tp=%d\n\tq=%d\n\tConstant factor: %f\n\n',muInd,etaInd,p,q,constFactor);
    end
    
    integral = 0;
    
    % Loop over vertices
    for i=1:nVerts
        
        % Get start and end indices for each edge
        indStart = rem(i-1,nVerts)+1;
        indEnd = rem(i,nVerts)+1;
        
        
%         [vertsProj(indStart,muInd), vertsProj(indEnd,muInd)]
%         [vertsProj(indStart,etaInd), vertsProj(indEnd,etaInd)]
%         [indStart, indEnd]
        
        % Evaluate contour integral along this line segment
        % evaluateSimpleLineIntegral(muStart,muEnd,etaStart,etaEnd,p,q)
        tmp = evaluateSimpleLineIntegral(vertsProj(indStart,muInd),  vertsProj(indEnd,muInd),  ...
                                         vertsProj(indStart,etaInd), vertsProj(indEnd,etaInd), ...
                                         vertsProj(indStart,differentialInd), vertsProj(indEnd,differentialInd), ...
                                         p, q);
        integral = integral + tmp;
                                                                    
    end
    
    % Unproject resulting integral
    % http://softsurfer.com/Archive/algorithm_0101/algorithm_0101.htm#3D Project to 2D
    areaScalingFactor = 1/(planeNormal(projectionPlane)/norm(planeNormal));
    res = constFactor*integral*areaScalingFactor;
    return;
    
    
    
    
    

% Do cylindrical stuff
elseif(strcmp(coordSys,'cylindrical'))
    
    
    % BetaGamma plane
    if(projectionPlane==1)
        muInd = 1;
        etaInd = 2;
        differentialInd = 3;
        p=1;
        q=1;
        constFactor = 1;
        fprintf('Choosing BG plane\n\tmuInd=%d\n\tetaInd=%d\n\tp=%d\n\tq=%d\n\tConstant factor: %f\n\n',muInd,etaInd,p,q,constFactor);
    % AlphaGamma plane
    elseif(projectionPlane==2)
        muInd = 3;
        etaInd = 1;
        differentialInd = 1;
        p=1;
        q=0;
        constFactor = 1;
        fprintf('Choosing AG plane\n\tmuInd=%d\n\tetaInd=%d\n\tp=%d\n\tq=%d\n\tConstant factor: %f\n\n',muInd,etaInd,p,q,constFactor);
    % AlphaBeta plane
    else
        muInd = 1;
        etaInd = 2;
        differentialInd = 2;
        p=2;
        q=0;
        constFactor = 0.5;
        fprintf('Choosing AB plane\n\tmuInd=%d\n\tetaInd=%d\n\tp=%d\n\tq=%d\n\tConstant factor: %f\n\n',muInd,etaInd,p,q,constFactor);
    end
    
    integral = 0;
    
    % Loop over vertices
    for i=1:nVerts
        
        % Get start and end indices for each edge
        indStart = rem(i-1,nVerts)+1;
        indEnd = rem(i,nVerts)+1;
        
        
%         [vertsProj(indStart,muInd), vertsProj(indEnd,muInd)]
%         [vertsProj(indStart,etaInd), vertsProj(indEnd,etaInd)]
%         [indStart, indEnd]
        
        % Evaluate contour integral along this line segment
        % evaluateSimpleLineIntegral(muStart,muEnd,etaStart,etaEnd,p,q)
        tmp = evaluateSimpleLineIntegral(vertsProj(indStart,muInd),  vertsProj(indEnd,muInd),  ...
                                         vertsProj(indStart,etaInd), vertsProj(indEnd,etaInd), ...
                                         vertsProj(indStart,differentialInd), vertsProj(indEnd,differentialInd), ...
                                         p, q);
        integral = integral + tmp;
                                                                    
    end
    
    % Unproject resulting integral
    % http://softsurfer.com/Archive/algorithm_0101/algorithm_0101.htm#3D Project to 2D
    areaScalingFactor = 1/(planeNormal(projectionPlane)/norm(planeNormal));
    res = constFactor*integral*areaScalingFactor;
    return;
    
    
    
% Do spherical stuff
elseif(strcmp(coordSys,'spherical'))
    
    % BetaGamma plane
    if(projectionPlane==1)
        muInd = 2;
        etaInd = 1;
        xiInd = 3;
        differentialInd = 3;
        p=1;
        q=0;
        constFactor = minComponent;%^2;
        
%         muInd = 2;
%         etaInd = 1;
%         xiInd = 3;
%         differentialInd = 2;
%         p=0;
%         q=0;
%         constFactor = minComponent^2;
        
        fprintf('Choosing BG plane\n\tmuInd=%d\n\tetaInd=%d\n\tp=%d\n\tq=%d\n\tConstant factor: %f\n\n',muInd,etaInd,p,q,constFactor);
    % AlphaGamma plane
    elseif(projectionPlane==2)
        muInd = 1;
        etaInd = 3;
        xiInd = 2;
        differentialInd = 3;
        p=2;
        q=0;
        constFactor = -(1/2);
        
        muInd = 1;
        etaInd = 3;
        xiInd = 2;
        differentialInd = 1;
        p=1;
        q=1;
        constFactor = 1;
        
        fprintf('Choosing AG plane\n\tmuInd=%d\n\tetaInd=%d\n\tp=%d\n\tq=%d\n\tConstant factor: %f\n\n',muInd,etaInd,p,q,constFactor);
    % AlphaBeta plane
    else
        muInd = 2;
        etaInd = 1;
        xiInd = 3;
        differentialInd = 1;
        p=1;
        q=1;
        constFactor = -sin(minComponent);%(1/2)*sin(minComponent);
        fprintf('Choosing AB plane\n\tmuInd=%d\n\tetaInd=%d\n\tp=%d\n\tq=%d\n\tConstant factor: %f\n\n',muInd,etaInd,p,q,constFactor);
    end
    
    integral = 0;
    
    % Loop over vertices
    for i=1:nVerts
        
        % Get start and end indices for each edge
        indStart = rem(i-1,nVerts)+1;
        indEnd = rem(i,nVerts)+1;
        
%         fprintf('-------- Edge %d ---------\n',i);
                
        % Evaluate contour integral along this line segment
        % If we're in the BetaGamma plane we need to call the complicated
        % integral
        % evaluateComplicatedLineIntegralSin(muStart,muEnd,etaStart,etaEnd,xiStart, xiEnd, differentialStart, differentialEnd, p,q)
        if(projectionPlane==1)
            tmp = evaluateComplicatedLineIntegralSin(vertsProj(indStart,muInd),  vertsProj(indEnd,muInd), ...
                                                     vertsProj(indStart,etaInd), vertsProj(indEnd,etaInd), ...
                                                     vertsProj(indStart,xiInd), vertsProj(indEnd,xiInd), ...
                                                     vertsProj(indStart,differentialInd), vertsProj(indEnd,differentialInd), ...
                                                     p,q);
        else
        % evaluateSimpleLineIntegral(muStart,muEnd,etaStart,etaEnd,p,q)
        tmp = evaluateSimpleLineIntegral(vertsProj(indStart,muInd),  vertsProj(indEnd,muInd),  ...
                                         vertsProj(indStart,etaInd), vertsProj(indEnd,etaInd), ...
                                         vertsProj(indStart,differentialInd), vertsProj(indEnd,differentialInd), ...
                                         p, q);
        end
        
        integral = integral + tmp;
                                                                    
    end
    
    % Unproject resulting integral
    % http://softsurfer.com/Archive/algorithm_0101/algorithm_0101.htm#3D Project to 2D
    areaScalingFactor = 1/(planeNormal(projectionPlane)/norm(planeNormal));
    res = constFactor*integral*areaScalingFactor;
    return;
    
    
    
% Error because the coordinate system is not supported
else
    fprintf(2,'[integrateContour] Coordinate system \"%s\" not supported\n',coordSys);
    return;
end