x = ones(3,1);
p = ones(3,1);
x = [1 0 0]
y = [0 1 0]
z = [0 0 1]
A = [x;y;z]
q = [1,1,1]
r = cross(q,rand(1,3));
s = cross(q,r);
q = q./(norm(q,2));
r = r./(norm(r,2));
s = s./(norm(s,2));
B = [q;r;s];
ex = cross(s,z)
ex = ex / (norm(ex,2))
theta = acos(dot(s,z))
ux  = ex(1);
uy = ex(2);
uz = ex(3);
C = cos(theta);
S = sin(theta);
t = 1 - C;
T = [t*ux^2+C t*ux*uy-S*uz t*ux*uz+S*uy; t*ux*uy+S*uz t*uy^2+C t*uy*uz-S*ux; t*ux*uz-S*uy t*uy*uz+S*ux t*uz^2+C];
rot1 = T*s';
rot2 = T*r';
rot3 = T*q';
Rtheta = [rot1,rot2,rot3]'

quiver3(p,p,p,A(:,1),A(:,2),A(:,3))
hold on
quiver3(p,p,p,Rtheta(:,1),Rtheta(:,2),Rtheta(:,3),'r')
quiver3(p,p,p,B(:,1),B(:,2),B(:,3),'g')
