clear all; close all; clc;

bnds = [0, 1, 0, 1, 0, 1];
nFaces = 6;

bbToVerts = [1,1,1;
             2,1,1;
             2,1,2;
             1,1,2;
             1,2,1;
             2,2,1;
             2,2,2;
             1,2,2];

verts = zeros(size(bbToVerts));
for i=1:8
    verts(i,:) = bnds(bbToVerts(i,:)+[0,2,4]);
end

quadFaces = [1,2,3,4;
             1,5,6,2;
             2,6,7,3;
             3,7,8,4;
             1,4,8,5;
             5,8,7,6];

polygons = cell(nFaces,1);

for i=1:nFaces
    
    polygons{i} = newPolygon();
    polygons{i}.vertices = quadFaces(i,:);
   
end
             
nPolygons = numel(polygons);

% polys = cell(1,1);
% polys(1) = polygons{1};

[volume, surfaceArea] = computeVolumeIntegrals(polygons(1:6), verts)