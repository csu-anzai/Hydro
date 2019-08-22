clear all; close all; clc;

nPts = 5;

         
planeNormal = [0,1,1]; 
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
vertsProj = projectPointsToPlane(verts,planeNormal);

area = 0;
for i=1:nPts
    
    indStart = rem(i-1,nPts)+1;
    indEnd = rem(i,nPts)+1;
    vStart = vertsProj(indStart,:);
    vEnd   = vertsProj(indEnd,:);
    dx = vEnd(1)-vStart(1);
    area = area + vStart(1)*vEnd(2) - vStart(2)*vEnd(1);

%     area = area + dx*(vStart(1)+0.5*dx);
    
    
end

area = area/2;

trueArea = 0.5*nPts*planeRadius^2*sin(2*pi/nPts);

(area-trueArea)/trueArea

% vertsABG = verts;
% vertsABG(:,3) = 0;
[T1, F1] = computeVolumeIntegrals(verts, planeNormal)


    
figure
plot(vertsProj(:,1),vertsProj(:,2))
xlim([-1,1])
ylim([-1,1])



% patch(planeXnew,planeYnew,planeZnew,'g','facealpha',0.05)


    
    