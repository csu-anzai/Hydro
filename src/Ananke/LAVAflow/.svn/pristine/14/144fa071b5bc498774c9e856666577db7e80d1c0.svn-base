function isInFront = verticesInfrontOfPlane(vertices, planeNormal, planePoint)

[nVerts, nComponents] = size(vertices);

isInFront = false(nVerts,1);

for i=1:nVerts
    vertices(i,:)
    diffPt = vertices(i,:)-planePoint;
    diffPt = diffPt/norm(diffPt);
    distanceFromPlane = dot(planeNormal,diffPt);
    
    if(distanceFromPlane>=0)
        isInFront(i) = true;
    end
end

	