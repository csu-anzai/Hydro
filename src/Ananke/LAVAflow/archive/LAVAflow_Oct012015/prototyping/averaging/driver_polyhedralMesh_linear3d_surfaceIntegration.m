clear all; close all; clc;

Xbnds = [0,1];
Ybnds = [0,2*pi];
Zbnds = [0,1];

Xbnds = [0,1];
Ybnds = [0,1];
Zbnds = [0,1];

% Xbnds = [0,1];
% Ybnds = [0,2*pi];
% Zbnds = [0,pi];

nCellX = 2;
nCellY = 3;
nCellZ = 2;
nCellTotal = nCellX*nCellY*nCellZ;

xBasis = linspace(Xbnds(1),Xbnds(2),nCellX+1);
yBasis = linspace(Ybnds(1),Ybnds(2),nCellY+1);
zBasis = linspace(Zbnds(1),Zbnds(2),nCellZ+1);

csCurrent = 'cartesian';

planeNormal = [1,0,0];
planePoint  = [0.5,0,0];


figure
hold on;
axis equal
area = 0;
numberOfIntegrations = 0;
for i=1:nCellX
    for j=1:nCellY
        for k=1:nCellZ
    
            XbndsCell = [xBasis(i),xBasis(i+1)];
            YbndsCell = [yBasis(j),yBasis(j+1)];
            ZbndsCell = [zBasis(k),zBasis(k+1)];
            [F,E,V] = create_box(XbndsCell,YbndsCell,ZbndsCell);
            
            
            vEdgeStart = V(E(:,1),:);
            vEdgeEnd  = V(E(:,2),:);

            pIntersections = intersectPlaneWithPolyhedra(vEdgeStart,vEdgeEnd, planePoint, planeNormal);
            nIntersections = size(pIntersections,1);
            
            if(nIntersections>2)
                tmp = integrateMonomialOverContour(csCurrent, pIntersections, planeNormal,1,0,0,0);
                area = area + tmp;
                numberOfIntegrations = numberOfIntegrations + 1;
            end
            
            
            % Drawing stuff from here down
            if(nIntersections>0)
                if(strcmp(csCurrent,'cartesian'))
                    % do nothing
                elseif(strcmp(csCurrent,'cylindrical'))
                    
                    pIntersections = [pIntersections(:,1).*cos(pIntersections(:,2)), ...
                                      pIntersections(:,1).*sin(pIntersections(:,2)), ...
                                      pIntersections(:,3)];
                else
                    pIntersections = [pIntersections(:,1).*cos(pIntersections(:,2)).*sin(pIntersections(:,3)), ...
                                      pIntersections(:,1).*sin(pIntersections(:,2)).*sin(pIntersections(:,3)), ...
                                      pIntersections(:,1).*cos(pIntersections(:,3))];
%                     plot3(ax,V(:,1).*cos(V(:,2)).*sin(V(:,3)),V(:,1).*sin(V(:,2)).*sin(V(:,3)), V(:,1).*cos(V(:,3)),'bo');
                end
                
                plot3(pIntersections(:,1),pIntersections(:,2),pIntersections(:,3),'rx')
            end
            if(nIntersections>2)
                nIntersections
            for ind=1:nIntersections
                if(nIntersections>2)
                    indStart = rem(ind-1,nIntersections)+1;
                    indEnd   = rem(ind,nIntersections)+1;
                end
                vStart   = pIntersections(indStart,:);
                vEnd     = pIntersections(indEnd,:);
%                 % transform vertices
%                 if(strcmp(csCurrent,'cartesian'))
%                     % do nothing
%                 elseif(strcmp(csCurrent,'cylindrical'))
%                     vStart = [vStart(1)*cos(vStart(2)),vStart(1)*sin(vStart(2)),vStart(3)];
%                     vEnd   = [vEnd(1)*cos(vEnd(2)),vEnd(1)*sin(vEnd(2)),vEnd(3)];
%                 end
                e = vEnd-vStart;
                quiver3(vStart(1),vStart(2),vStart(3),e(1),e(2),e(3),'g','linewidth',2);
            end
            end
            
            draw_box(gca,F,E,V,csCurrent)
        end
    end
    
end

fprintf('area = %f\n',area);

%%
% if(strcmp(csCurrent,'cartesian'))
% % do nothing
% elseif(strcmp(csCurrent,'cylindrical'))
% 
%     planePoint = [planePoint(1)*cos(planePoint(2)),planePoint(1)*sin(planePoint(2)),planePoint(3)];
%     planeNormal = [planeNormal(1)*cos(planePoint(2)),planeNormal(2)*sin(planePoint(2)),planeNormal(3)];
%     
% end


if(strcmp(csCurrent,'cartesian'))
    % do nothing
elseif(strcmp(csCurrent,'cylindrical'))
    planeNormal = [planeNormal(1)*cos(planePoint(2)),planeNormal(2)*sin(planePoint(2)),planeNormal(3)];
end

% planeTheta = linspace(0,2*pi,5); planeTheta = planeTheta(1:end-1);
% planeRadius = 1.5;
% planeX = planeRadius*cos(planeTheta);
% planeY = planeRadius*sin(planeTheta);
% planeZ = zeros(size(planeX));
% tmpBasis = zeros(3,1); [~,ind] = min(abs(planeNormal));tmpBasis(ind)=1;
% uBasis = cross(planeNormal,tmpBasis); uBasis = uBasis/norm(uBasis);
% vBasis = cross(planeNormal,uBasis);   vBasis = vBasis/norm(vBasis);
% planeXnew = planeX*uBasis(1) + planeY*vBasis(1) + planePoint(1);
% planeYnew = planeX*uBasis(2) + planeY*vBasis(2) + planePoint(2);
% planeZnew = planeX*uBasis(3) + planeY*vBasis(3) + planePoint(3);
% 
% planeSamples = [planeXnew; planeYnew; planeZnew]';
% 
% 
% patch(planeSamples(:,1),planeSamples(:,2),planeSamples(:,3),'g','facealpha',0.5)
