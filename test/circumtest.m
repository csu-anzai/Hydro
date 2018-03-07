clc, clear all, clf

% define triangle vertices
x = [0 .5 1];
y = [0 2 .2];
cor = [x;y];

v1 = [0 0];
v2 = [.5 2];
v3 = [1 .2];

[x y] = circumcenter(v1,v2,v3)

% compute plot
[r,cn] = circumcircle(cor,1)
