clc, clear all

% number of grid points in each direction
numPts = 6;

% error threshold
threshold = 1e-2;

% bounding box
a = 0.0;
b = 1;

% create linear arrays
x = linspace(a,b,numPts);
y = linspace(a,b,numPts);

% create grid array
[X,Y] = meshgrid(x,y);

% add perturbation to interior points
eps = 1e-1;
for i=2:length(X)-1
    for j = 2:length(Y)-1
        X(i,j) = X(i,j) + 2*eps*rand - eps;
        Y(i,j) = Y(i,j) + 2*eps*rand - eps;
    end
end

% create suitable array for matlab built-in function
k = 1;
for i=1:length(X)
    for j=1:length(X)
        P(k,1) = X(i,j);
        P(k,2) = Y(i,j);
        k = k + 1;
    end
end

% compute the delaunay triangulation
DT = delaunayTriangulation(P)

% define functional data at points
f = sin(P(:,1)) .* cos(P(:,2));

% evaluate odd points one-by-one for potential removal

% compute the approximate nodes approximate functional value via natural neighbors
vq = griddata(P(:,1),P(:,2),f,xq,yq)

% compute actual value based on known function
vtrue = sin(xq) * cos(yq)

% compute the detail coefficient
d = vtrue - vq

% plot the triangulation
figure(1)
triplot(DT);
