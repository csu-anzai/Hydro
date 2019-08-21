function vertsProj = projectPointsToABGPlane(verts,polyNormal,origin)


nVerts = size(verts,1);

% Get the right-handed coordinate system that maximizes the 3rd component
% of the normal in the new system
csPerm = getCoordSystemPermutation(polyNormal);

% Permute the points to the proper new coordinate system
vertsPerm = verts(:,csPerm);
polyNormalPerm = polyNormal(csPerm);
ABPlane = polyNormal;
vertsProj = zeros(size(verts));
for i=1:nVerts
%     tmp = origin - vertsPerm(i,:);
    tmp = origin - verts(i,:);
    tmp = tmp - dot(ABPlane,tmp)*ABPlane;
    vertsProj(i,:) = origin-tmp;
end


