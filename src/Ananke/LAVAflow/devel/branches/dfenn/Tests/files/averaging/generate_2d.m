clear all; close all; clc;

filename = 'structpts_constant_2d.vtk';
nCellX = 3; nCellY = 3; nCellZ = 1;

xbnds = [0, 1];
ybnds = [-2,2];
zbnds = zeros(1,2);

xCellFaces = linspace(xbnds(1),xbnds(2),nCellX+1);
yCellFaces = linspace(ybnds(1),ybnds(2),nCellY+1);

dx = diff(xCellFaces); dx = dx(1);
dy = diff(yCellFaces); dy = dy(1);
dz = 0;

xCellCenters = linspace(xbnds(1)+dx/2, xbnds(2)-dx/2, nCellX);
yCellCenters = linspace(ybnds(1)+dy/2, ybnds(2)-dy/2, nCellY);

[X,Y] = meshgrid(xCellCenters, yCellCenters);
             
F1 = ones(size(X)); % Constant field
F2 = X; % Constant field



fid = fopen(filename,'w');

fprintf(fid,'# vtk DataFile Version 2.0\n');
fprintf(fid,'Structured points constant data\n');
fprintf(fid,'ASCII\n');
fprintf(fid,'DATASET RECTILINEAR_GRID\n');
fprintf(fid,'DIMENSIONS %d %d %d\n',nCellX+1, nCellY+1, nCellZ+1);
fprintf(fid,'X_COORDINATES %d float\n',nCellX+1);
for i=1:numel(xCellFaces);
    fprintf(fid,'%f\n',xCellFaces(i));
end

fprintf(fid,'Y_COORDINATES %d float\n',nCellY+1);
for i=1:numel(yCellFaces);
    fprintf(fid,'%f\n',yCellFaces(i));
end

fprintf(fid,'Z_COORDINATES %d float\n%f\n%f\n',nCellZ+1,0,0);

fprintf(fid, '\n');
fprintf(fid, 'CELL_DATA   %d\n', nCellX*nCellY);


%fprintf(fid, 'CELL_DATA   %d\n', nCellX*nCellY);
fprintf(fid, 'SCALARS vari2 float\n');
fprintf(fid, 'LOOKUP_TABLE default\n');
for a=1:nCellX
    for b=1:nCellY
        fprintf(fid, '%f\n', F2(a,b));
    end
%     fprintf(fid, '\n');
end


fprintf(fid, '\n');
fprintf(fid, 'SCALARS vari1 float\n');
fprintf(fid, 'LOOKUP_TABLE default\n');
for a=1:nCellX
    for b=1:nCellY
        fprintf(fid, '%f\n', F1(a,b));
    end
%     fprintf(fid, '\n');
end


fclose(fid);