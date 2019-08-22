clear all; close all; clc;

%% Test AB plane (sides of cones)
% Test AB plane (sides of cones)
planeNormal = [0,0,1];
rInner = 0.5;
rOuter = 1.0;
thetaMin = 0;
thetaMax = 2*pi;
phi = pi/3;

verts = [rInner, thetaMin, phi;
         rOuter, thetaMin, phi;
         rOuter, thetaMax, phi;
         rInner, thetaMax, phi];

nVerts = size(verts,1);


trueArea = 0.5*(thetaMax-thetaMin)*sin(phi)*(rOuter^2-rInner^2);

% area = integrateContour('spherical', verts, planeNormal);
area = integrateMonomialOverContour('spherical', verts, planeNormal,1,0,0,0);


fprintf('True Area: %5.5E\nComputed Area: %5.5E\nRel Error: %5.5E\n',trueArea,area,(area-trueArea)/trueArea);


    
figure
ptsPerLine = 100;
vertsResamp = [];
for i=1:nVerts
    
    indStart = rem(i-1,nVerts)+1;
    indEnd = rem(i,nVerts)+1;
        
    r = linspace(verts(indStart,1),verts(indEnd,1),ptsPerLine);
    theta = linspace(verts(indStart,2),verts(indEnd,2),ptsPerLine);
    %phi = linspace(verts(indStart,3),verts(indEnd,3),ptsPerLine);
    vertsResamp = [vertsResamp;r',theta'];
end

[R,T] = meshgrid(vertsResamp(:,1),vertsResamp(:,2));

x = R.*cos(T).*sin(phi);
y = R.*sin(T).*sin(phi);
z = R.*cos(phi);
mesh(x,y,z);
% patch(x,y,z,'r')
% plot(x,y)
axis equal
xlim([-1,1])
ylim([-1,1])
zlim([-1,1])


%% Test AG plane (disks formed by great circles)
% Test AG plane (disks formed by great circles)
planeNormal = [0,1,0];
rInner = 0.5;
rOuter = 1.0;
theta = sqrt(2);
phiMin = 0;
phiMax = pi;

verts = [rInner, theta, phiMin;
         rInner, theta, phiMax;
         rOuter, theta, phiMax;
         rOuter, theta, phiMin];

nVerts = size(verts,1);


trueArea = (phiMax-phiMin)/2*(rOuter^2-rInner^2);

% area = integrateContour('spherical', verts, planeNormal);
area = integrateMonomialOverContour('spherical', verts, planeNormal,1,0,0,0);


fprintf('True Area: %5.5E\nComputed Area: %5.5E\nRel Error: %5.5E\n',trueArea,area,(area-trueArea)/trueArea);


    
figure
ptsPerLine = 100;
vertsResamp = [];
for i=1:nVerts
    
    indStart = rem(i-1,nVerts)+1;
    indEnd = rem(i,nVerts)+1;
        
    r = linspace(verts(indStart,1),verts(indEnd,1),ptsPerLine);
    phi = linspace(verts(indStart,3),verts(indEnd,3),ptsPerLine);
    %phi = linspace(verts(indStart,3),verts(indEnd,3),ptsPerLine);
    vertsResamp = [vertsResamp;r',phi'];
end

[R,P] = meshgrid(vertsResamp(:,1),vertsResamp(:,2));

x = R.*cos(theta).*sin(P);
y = R.*sin(theta).*sin(P);
z = R.*cos(P);
mesh(x,y,z);
% patch(x,y,z,'r')
% plot(x,y)
axis equal
xlim([-1,1])
ylim([-1,1])
zlim([-1,1])

%% Test BG plane (patch on the surface of the sphere)
% Test BG plane (patch on the surface of the sphere)
planeNormal = [1,0,0];
R = 1;
thetaMin = 0;
thetaMax = 0.5*pi;
phiMin = .4*pi;
phiMax = .6*pi;

verts = [R, thetaMin, phiMin;
         R, thetaMax, phiMin;
         R, thetaMax, phiMax;
         R, thetaMin, phiMax];

nVerts = size(verts,1);


trueArea = R^2*(thetaMax-thetaMin)*(cos(phiMin)-cos(phiMax));

% area = integrateContour('spherical', verts, planeNormal);
area = integrateMonomialOverContour('spherical', verts, planeNormal,1,0,0,0);


fprintf('True Area: %5.5E\nComputed Area: %5.5E\nRel Error: %5.5E\n',trueArea,area,(area-trueArea)/trueArea);


    
figure
ptsPerLine = 100;
vertsResamp = [];
for i=1:nVerts
    
    indStart = rem(i-1,nVerts)+1;
    indEnd = rem(i,nVerts)+1;
        
    theta = linspace(verts(indStart,2),verts(indEnd,2),ptsPerLine);
    phi = linspace(verts(indStart,3),verts(indEnd,3),ptsPerLine);
    %phi = linspace(verts(indStart,3),verts(indEnd,3),ptsPerLine);
    vertsResamp = [vertsResamp;theta',phi'];
end

[T,P] = meshgrid(vertsResamp(:,1),vertsResamp(:,2));

x = R.*cos(T).*sin(P);
y = R.*sin(T).*sin(P);
z = R.*cos(P);
mesh(x,y,z);
% patch(x,y,z,'r')
% plot(x,y)
axis equal
xlim([-1,1])
ylim([-1,1])
zlim([-1,1])

