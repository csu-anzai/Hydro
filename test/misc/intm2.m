function r = intm2(f, t)
% f: function 
% t: three points of a triangle
% r: integration of f over t
a = t(1,:);
b = t(2,:);
c = t(3,:);

jg = abs((b(1)-a(1))*(c(2)-a(2))-(c(1)-a(1))*(b(2)-a(2)));
ftilda = @(u,v) f(a(1)+u*(b(1)-a(1))+v*(c(1)-a(1)), a(2)+u*(b(2)-a(2))+v*(c(2)-a(2)));

fy = @(x) 1-x;
r = jg * integral2(ftilda, 0,1, 0,fy,'AbsTol',1e-15);

end
