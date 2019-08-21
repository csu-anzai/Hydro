function normals = getPolygonNormals(polygons, verts)

nPolys = length(polygons);
normals = zeros(nPolys,3);

for i=1:nPolys
    e1 = verts(polygons{i}.vertices(2),:)-verts(polygons{i}.vertices(1),:);
    e2 = verts(polygons{i}.vertices(3),:)-verts(polygons{i}.vertices(2),:);
    normals(i,:) = cross(e1,e2);
end