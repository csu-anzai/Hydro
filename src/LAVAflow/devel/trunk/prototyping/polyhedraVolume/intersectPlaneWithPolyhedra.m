function pIntersections = intersectPlaneWithPolyhedra(vEdgeStart,vEdgeEnd, pPlane, normalPlane)

[~,maxInd] = max(abs(normalPlane));
projInds = [1:maxInd-1,maxInd+1:3];


% Loop over edges
% vEdgeStart = V(E(:,1),:);
% vEdgeEnd  = V(E(:,2),:);
pIntersections = [];
for e=1:size(vEdgeStart,1)
    pStart = vEdgeStart(e,:);
    pEnd = vEdgeEnd(e,:);
    [isIntersecting, ~, pInter] = linePlaneIntersection(pStart, pEnd, pPlane, normalPlane);
    if(isIntersecting)
        pIntersections(end+1,:) = pInter;
    end
end

if(numel(pIntersections)==0)
    return;
end


% Sometimes the plane may cut through a corner. This will, essentially,
% duplicate vertices
% Clean this problem up now
pIntersectionsNew = pIntersections(1,:);
for i=1:size(pIntersections,1)
    inB = false;
    for j=1:size(pIntersectionsNew,1)
        if(norm(pIntersections(i,:)-pIntersectionsNew(j,:))<1e-6)
            inB = true;
            break;
        end
    end
    if(~inB)
        pIntersectionsNew(end+1,:) = pIntersections(i,:);
    end
end

pIntersections = pIntersectionsNew;
nIntersections = size(pIntersections,1);
            


% Project the points into a 2d plane
if(nIntersections>2)
    pProj = pIntersections(:,projInds);
    pProjCenter = mean(pProj);
    pProj = bsxfun(@minus,pProj,pProjCenter);
    % Determine the angle from the new x axis
    theta = atan2(pProj(:,2),pProj(:,1));
    % sort to get order
    [~,sortInd] = sort(theta);

    pIntersections = pIntersections(sortInd,:);
end

% The order may be incorrect, so compute the normal from the points and
% reverse the order if needed
if(nIntersections>2)
    computedNormal = cross(pIntersections(2,:)-pIntersections(1,:),pIntersections(3,:)-pIntersections(2,:));
    computedNormal = computedNormal/norm(computedNormal);
    if(dot(computedNormal,normalPlane)<0)
        pIntersections = flipud(pIntersections);
    end
end