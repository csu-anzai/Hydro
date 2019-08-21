clear all; close all; clc;

verts = [1, 0;
         1,pi/2;
         0,pi/2;
         0,0];
     
nPts = size(verts,1);
    
area = 0;

for i=1:nPts
    
    indStart = rem(i-1,nPts)+1;
    indEnd = rem(i,nPts)+1;
    vStart = verts(indStart,:);
    vEnd   = verts(indEnd,:);
    
    dr = vEnd(1)-vStart(1);
    dT = vEnd(2)-vStart(2);
    r0T0 = vStart(1)*vStart(2);
    r1T1 = vEnd(1)*vEnd(2);
    
    area = area + -dr*((1/3)*dr*dT + (1/2)*(r1T1-r0T0-dr*dT) + r0T0);

    
end

area
areaTrue = pi/4*(1-0^2)


plot(verts(:,1).*cos(verts(:,2)),verts(:,1).*sin(verts(:,2)))