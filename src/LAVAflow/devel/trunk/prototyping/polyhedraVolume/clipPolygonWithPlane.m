function [clippedPoly, vertsNew,intersectionPoints] = clipPolygonWithPlane(origPoly, verts, planeNormal, planePoint)

isInFront = verticesInfrontOfPlane(verts,planeNormal,planePoint);

nIntersections = 0;
intersectionPoints = zeros(1,3);
foundFirstCut = false;
foundLastCut  = false;

nFaceVerts = length(origPoly.vertices);
nVerts = size(verts,1);

verticesNew = [];


% Vertex to start traversing over
ind0 = find(isInFront,1);

% Loop over all edges
for j=ind0:nFaceVerts+ind0-1
    indStart = rem(j-1,nFaceVerts)+1;
    indEnd = rem(j,nFaceVerts)+1;
    vStart = verts(origPoly.vertices(indStart),:);
    vEnd   = verts(origPoly.vertices(indEnd),:);

    t = dot(planeNormal, planePoint-vStart)/dot(planeNormal,vEnd-vStart);

    if(~foundFirstCut && ~foundLastCut)
        verticesNew(end+1) = origPoly.vertices(indStart);
    end

    % Edge intersects plane
    if(t > 0 && t< 1)
        nIntersections = nIntersections + 1;
        intersectionPoints(nIntersections,:) = vStart + t*(vEnd-vStart);
        
        if(~foundFirstCut)
            foundFirstCut = true;
            verticesNew(end+1) = nVerts + nIntersections;
        else
            foundLastCut = true;
            verticesNew(end+1) = nVerts + nIntersections;

        end
    end

end

clippedPoly = newPolygon();
clippedPoly.vertices = verticesNew;    

vertsNew = [verts;intersectionPoints];