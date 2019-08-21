clear all; close all; clc;

bnds = [0, 1, 0, 1, 0, 1];
nFaces = 1;

bbToVerts = [1,1,1;
             2,1,1;
             2,2,1;
             1,2,1];

verts = zeros(size(bbToVerts));
nVerts = 4;
for i=1:nVerts
    verts(i,:) = bnds(bbToVerts(i,:)+[0,2,4]);
end

quadFaces = [1,2,3,4];
             

polygons = cell(nFaces,1);

for i=1:nFaces
    
    polygons{i} = newPolygon();
    polygons{i}.vertices = quadFaces(i,:);
   
end
             
nPolygons = numel(polygons);
         
planeNormal = [1,1,0]; planeNormal = planeNormal/norm(planeNormal);
planePoint  = [0.25,0.5,0];

isInFront = verticesInfrontOfPlane(verts,planeNormal,planePoint);

[polygons{1}, verts, intersectionPoints] = clipPolygonWithPlane(polygons{1}, verts, planeNormal, planePoint)


polyNormals = getPolygonNormals(polygons,verts);



figure
hold on;
for i=1:nPolygons
    
    patch(verts(polygons{i}.vertices,1), ...
          verts(polygons{i}.vertices,2), ...
          verts(polygons{i}.vertices,3), ...
          'r','facealpha',0.1,'linewidth',1);
    
    nFaceVerts = length(polygons{i}.vertices);
    
    for j=1:nFaceVerts
        indStart = rem(j-1,nFaceVerts)+1;
        indEnd = rem(j,nFaceVerts)+1;
        vStart = verts(polygons{i}.vertices(indStart),:);
        vEnd   = verts(polygons{i}.vertices(indEnd),:);
        e = vEnd-vStart;
        quiver3(vStart(1),vStart(2),vStart(3),e(1),e(2),e(3),'g','linewidth',2);
    end
      
end


plot3(verts(isInFront,1),verts(isInFront,2),verts(isInFront,3),'ro')
quiver3(planePoint(1),planePoint(2),planePoint(3),planeNormal(1),planeNormal(2),planeNormal(3))

planeTheta = linspace(0,2*pi,50); planeTheta = planeTheta(1:end-1);
planeRadius = 1;
planeX = planeRadius*cos(planeTheta);
planeY = planeRadius*sin(planeTheta);
planeZ = zeros(size(planeX));

tmpBasis = zeros(3,1); [~,ind] = min(abs(planeNormal));tmpBasis(ind)=1;
uBasis = cross(planeNormal,tmpBasis); uBasis = uBasis/norm(uBasis);
vBasis = cross(planeNormal,uBasis);   vBasis = vBasis/norm(vBasis);
planeXnew = planeX*uBasis(1) + planeY*vBasis(1) + planePoint(1);
planeYnew = planeX*uBasis(2) + planeY*vBasis(2) + planePoint(2);
planeZnew = planeX*uBasis(3) + planeY*vBasis(3) + planePoint(3);

patch(planeXnew,planeYnew,planeZnew,'g','facealpha',0.05)

plot3(intersectionPoints(:,1),intersectionPoints(:,2),intersectionPoints(:,3),'rx')

axis equal

xlim([-1,2]);
ylim([-1,2]);
zlim([-1,2]);