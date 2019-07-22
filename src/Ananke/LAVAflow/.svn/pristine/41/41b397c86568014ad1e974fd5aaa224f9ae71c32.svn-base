clear all; close all; clc;

%% Test AB plane (r-theta disks)
% Test AB plane (r-theta disks)
planeNormal = [0,0,1];
rInner = 0.5;
rOuter = 1.0;
thetaMin = 0;
thetaMax = 2*pi;

verts = [rInner, thetaMin, 0;
         rOuter, thetaMin, 0;
         rOuter, thetaMax, 0;
         rInner, thetaMax, 0];

nVerts = size(verts,1);


trueArea = ((thetaMax-thetaMin)/2)*(rOuter^2-rInner^2);

% area = integrateContour('cylindrical', verts, planeNormal);
area = integrateMonomialOverContour('cylindrical', verts, planeNormal,1,0,0,0);

fprintf('True Area: %5.5E\nComputed Area: %5.5E\nRel Error: %5.5E\n',trueArea,area,(area-trueArea)/trueArea);


    
figure
ptsPerLine = 100;
vertsResamp = [];
for i=1:nVerts
    
    indStart = rem(i-1,nVerts)+1;
    indEnd = rem(i,nVerts)+1;
        
    r = linspace(verts(indStart,1),verts(indEnd,1),ptsPerLine);
    theta = linspace(verts(indStart,2),verts(indEnd,2),ptsPerLine);
    vertsResamp = [vertsResamp;r',theta'];
end

x = vertsResamp(:,1).*cos(vertsResamp(:,2));
y = vertsResamp(:,1).*sin(vertsResamp(:,2));
patch(x,y,'r')
% plot(x,y)
axis equal

%% Test AG plane (r-z rectangles)
planeNormal = [0,1,0];
rInner = 0.5;
rOuter = 1.0;
thetaMin = 0;
thetaMax = 0;
zMin = 0.0;
zMax = 2.0;

% CCW list
verts = [rInner, thetaMin, zMin;
         rInner, thetaMin, zMax;
         rOuter, thetaMax, zMax;
         rOuter, thetaMax, zMin];

     
nVerts = size(verts,1);


trueArea = (rOuter-rInner)*(zMax-zMin);

% area = integrateContour('cylindrical', verts, planeNormal);
area = integrateMonomialOverContour('cylindrical', verts, planeNormal,1,0,0,0);


fprintf('True Area: %5.5E\nComputed Area: %5.5E\nRel Error: %5.5E\n',trueArea,area,(area-trueArea)/trueArea);


%% Test BG plane (theta-z tubes)
planeNormal = [1,0,0];
R = 0.5;
thetaMin = 0;
thetaMax = 2*pi;
zMin = 0.0;
zMax = 2.0;

% CCW list
verts = [R, thetaMin, zMin;
         R, thetaMax, zMin;
         R, thetaMax, zMax;
         R, thetaMin, zMax];



trueArea = (thetaMax-thetaMin)*R*(zMax-zMin);

% area = integrateContour('cylindrical', verts, planeNormal);
area = integrateMonomialOverContour('cylindrical', verts, planeNormal,1,0,0,0);


fprintf('True Area: %5.5E\nComputed Area: %5.5E\nRel Error: %5.5E\n',trueArea,area,(area-trueArea)/trueArea);


figure
ptsPerLine = 100;
vertsResamp = [];
for i=1:nVerts
    
    indStart = rem(i-1,nVerts)+1;
    indEnd = rem(i,nVerts)+1;
        
    r = linspace(verts(indStart,1),verts(indEnd,1),ptsPerLine);
    theta = linspace(verts(indStart,2),verts(indEnd,2),ptsPerLine);
    z = linspace(verts(indStart,3),verts(indEnd,3),ptsPerLine);
    vertsResamp = [vertsResamp;r',theta',z'];
end

x = vertsResamp(:,1).*cos(vertsResamp(:,2));
y = vertsResamp(:,1).*sin(vertsResamp(:,2));
z = vertsResamp(:,3);
patch(x,y,z,'r')
% plot(x,y)
axis equal

