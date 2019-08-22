function [volume, surfaceArea] = computeVolumeIntegrals(polygons, verts)


nPolygons = numel(polygons);
polyNormals = getPolygonNormals(polygons, verts);

volume = 0;
surfaceArea = 0;

for i=1:nPolygons
    polyNormal = polyNormals(i,:);
    csPerm = getCoordSystemPermutation(polyNormal);
    
    polyNormalTrans = polyNormal(csPerm);
    polyVertsTrans = verts(polygons{i}.vertices,csPerm);

    [F1, Fa,Fb,Fg] = computeFaceIntegrals(polyVertsTrans, polyNormalTrans);
    
    surfaceArea = surfaceArea + abs(F1);
    
    if(csPerm(1)==1)
        volume = volume + polyNormalTrans(1)*Fa;
    elseif(csPerm(1)==2)
        volume = volume + polyNormalTrans(2)*Fb;
    else
        volume = volume + polyNormalTrans(3)*Fg;
    end
    
end