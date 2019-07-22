function verts = unprojectPointsFromPlane(vertsProj,planeNormal)

[~,ind] = max(abs(planeNormal));
minComp = min(vertsProj(:,ind));

verts = vertsProj;
verts(:,ind) = minComp;
