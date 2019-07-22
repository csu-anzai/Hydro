function vertsProj = projectPointsToPlane(verts,planeNormal)


nVerts = size(verts,1);


tmpBasis = zeros(3,1); [~,ind] = min(abs(planeNormal));tmpBasis(ind)=1;
uBasis = cross(planeNormal,tmpBasis); uBasis = uBasis/norm(uBasis);
vBasis = cross(planeNormal,uBasis);   vBasis = vBasis/norm(vBasis);
wBasis = planeNormal; wBasis = wBasis/norm(wBasis);

for i=1:nVerts
    vertsProj(i,1) = verts(i,1)*uBasis(1) + verts(i,2)*uBasis(2)  + verts(i,3)*uBasis(3);
    vertsProj(i,2) = verts(i,1)*vBasis(1) + verts(i,2)*vBasis(2)  + verts(i,3)*vBasis(3);
    vertsProj(i,3) = verts(i,1)*wBasis(1) + verts(i,2)*wBasis(2)  + verts(i,3)*wBasis(3);
end



