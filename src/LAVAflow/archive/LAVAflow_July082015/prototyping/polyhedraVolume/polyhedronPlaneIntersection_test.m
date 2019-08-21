clear all; close all; clc;

[F,E,V] = create_box([0,1],[0,1],[0,1]);
nFaces = size(F,1);
nEdges = size(E,1);

vEdgeStart = V(E(1:4,1),:);
vEdgeEnd  = V(E(1:4,2),:);

normalPlane = [1,1,1]; normalPlane = normalPlane/norm(normalPlane);
pPlane = [0.15,0.15,0.15];

pIntersections = intersectPlaneWithPolyhedra(vEdgeStart,vEdgeEnd, pPlane, normalPlane);
nIntersections = size(pIntersections,1);


if(nIntersections>2)
    computedNormal = cross(pIntersections(2,:)-pIntersections(1,:),pIntersections(3,:)-pIntersections(2,:));
    computedNormal = computedNormal/norm(computedNormal);
end

figure
axis equal;
modifyaxis
draw_box(gca,F,E,V,'cartesian')

hold on;

for j=1:nIntersections
    if(nIntersections>2)
        indStart = rem(j-1,nIntersections)+1;
        indEnd   = rem(j,nIntersections)+1;
    else
        indStart = j;
        indEnd = j+1;
        if(indEnd>nIntersections)
            break;
        end
    end
    vStart   = pIntersections(indStart,:);
    vEnd     = pIntersections(indEnd,:);
    e = vEnd-vStart;
    quiver3(vStart(1),vStart(2),vStart(3),e(1),e(2),e(3),'g','linewidth',2);
end

if(nIntersections>0)
    plot3(pIntersections(:,1),pIntersections(:,2),pIntersections(:,3),'rx')
end

planeTheta = linspace(0,2*pi,50); planeTheta = planeTheta(1:end-1);
planeRadius = 1;
planeX = planeRadius*cos(planeTheta);
planeY = planeRadius*sin(planeTheta);
planeZ = zeros(size(planeX));
tmpBasis = zeros(3,1); [~,ind] = min(abs(normalPlane));tmpBasis(ind)=1;
uBasis = cross(normalPlane,tmpBasis); uBasis = uBasis/norm(uBasis);
vBasis = cross(normalPlane,uBasis);   vBasis = vBasis/norm(vBasis);
planeXnew = planeX*uBasis(1) + planeY*vBasis(1) + pPlane(1);
planeYnew = planeX*uBasis(2) + planeY*vBasis(2) + pPlane(2);
planeZnew = planeX*uBasis(3) + planeY*vBasis(3) + pPlane(3);

patch(planeXnew,planeYnew,planeZnew,'g','facealpha',0.05)


quiver3(pPlane(1),pPlane(2),pPlane(3),normalPlane(1),normalPlane(2),normalPlane(3),'g','linewidth',2);
if(nIntersections>2)
    quiver3(pPlane(1),pPlane(2),pPlane(3),computedNormal(1),computedNormal(2),computedNormal(3),'b.','linewidth',2);
end
