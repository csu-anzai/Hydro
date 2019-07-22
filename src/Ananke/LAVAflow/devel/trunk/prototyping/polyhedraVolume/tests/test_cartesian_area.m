clear all; close all; clc;

nPts = 5;

         
planeNormal = [rand(1),rand(1),rand(1)]; 
planePoint  = [0,0,0];
planeNormal = planeNormal/norm(planeNormal);

planeTheta = linspace(0,2*pi,nPts+1); planeTheta = planeTheta(1:end-1);
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




verts = [planeXnew', planeYnew', planeZnew'];

trueArea = 0.5*nPts*planeRadius^2*sin(2*pi/nPts);

% area = integrateContour('cartesian', verts, planeNormal);
area = integrateMonomialOverContour('cartesian', verts, planeNormal,1,0,0,0);


fprintf('True Area: %5.5E\nComputed Area: %5.5E\nRel Error: %5.5E\n',trueArea,area,(area-trueArea)/trueArea);


    
figure
plot3([verts(:,1);verts(1,1)],[verts(:,2);verts(1,2)],[verts(:,3);verts(1,3)])
xlim([-1,1])
ylim([-1,1])