x = linspace(-2, 2, 200);
y = linspace(-2, 2, 200);
z = linspace(-2, 2, 200);
[ xx, yy, zz ] = meshgrid(x, y, z);
c = (xx).^2 + (yy).^2 + (zz).^2;
[T, p] = MarchingCubes(xx, yy, zz, c, 0.5);
trimesh(T, p(:, 1), p(:, 2), p(:, 3));curvature;
meancurv = mean(verteigs);



trisurf(T,p(:,1),p(:,2),p(:,3),meancurv')
colorbar