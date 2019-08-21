function draw_box(ax,F,E,V,coordSys)

doDrawVertices = true;
doDrawEdges    = true;
doDrawSides    = false;

CS_CART   = 1;
CS_CYL    = 2;
CS_SPHERE = 3;
csCurrent = 1;

if(strcmp(coordSys,'cartesian'))
    csCurrent = CS_CART;
elseif(strcmp(coordSys,'cylindrical'))
    csCurrent = CS_CYL;
elseif(strcmp(coordSys,'spherical'))
    csCurrent = CS_SPHERE;
else        
    fprintf(2,'[draw_box] Coordinate system \"%s\" not supported\n!',coordSys);
    return;
end




F=F;
E=E;

holdStatus = get(gca,'NextPlot');
set(gca,'NextPlot', 'add');

% Draw points
if(doDrawVertices)
    switch csCurrent
        case CS_CART
            plot3(ax,V(:,1),V(:,2),V(:,3),'bo');
        case CS_CYL
            plot3(ax,V(:,1).*cos(V(:,2)),V(:,1).*sin(V(:,2)), V(:,3),'bo');
        case CS_SPHERE
            plot3(ax,V(:,1).*cos(V(:,2)).*sin(V(:,3)),V(:,1).*sin(V(:,2)).*sin(V(:,3)), V(:,1).*cos(V(:,3)),'bo');
    end
end


% Draw edges
nPtsPerEdge = 20;
if(doDrawEdges)
    for i=1:size(E,1)
        pts = V(E(i,:),:);
        
        switch csCurrent
            case CS_CART
            
                line(pts(:,1),pts(:,2),pts(:,3))
%                 ray = pts(2,:)-pts(1,:);
%                 quiver3(pts(1,1),pts(1,2),pts(1,3), ray(1),ray(2),ray(3),'g','linewidth',2);
                
            
            case CS_CYL
                pStart = pts(1,:);
                pEnd   = pts(2,:);
                r     = linspace(pStart(1),pEnd(1),nPtsPerEdge);
                theta = linspace(pStart(2),pEnd(2),nPtsPerEdge);
                z     = linspace(pStart(3),pEnd(3),nPtsPerEdge);
                x = r.*cos(theta);
                y = r.*sin(theta);

                plot3(ax,x,y,z);

            case CS_SPHERE
                
                pStart = pts(1,:);
                pEnd   = pts(2,:);
                r     = linspace(pStart(1),pEnd(1),nPtsPerEdge);
                theta = linspace(pStart(2),pEnd(2),nPtsPerEdge);
                phi   = linspace(pStart(3),pEnd(3),nPtsPerEdge);
                
                x = r.*cos(theta).*sin(phi);
                y = r.*sin(theta).*sin(phi);
                z = r.*cos(phi);
                
                plot3(ax,x,y,z);
            otherwise
                fprintf(2,'[draw_box] Coordinate system \"%s\" not supported\n!',coordSys);
                return;
        end
        
        
    end
        
end


nPtsPerEdgeFace = 5;
if(doDrawSides)
    
    nFaces = size(F,1);
    
    
    switch csCurrent
        case CS_CART
            for f=1:nFaces
                edges = E(F(f,:),:);
%                 edges = unique(edges);
                nEdges = size(edges,1);
                verts = V(edges,:)
                
                planeNormal = cross(verts(4,:)-verts(2,:),verts(3,:)-verts(1,:));
                planeNormal = planeNormal/norm(planeNormal);
                
                [~,c] = max(abs(planeNormal));
                minComp = min(verts(:,c));
                newPlaneNormal = [0, 0, 0];
                newPlaneNormal(c) = sign(planeNormal(c));
                
                inds = [1:c-1, c+1:3];
%                 planeNormal
%                 inds
                vertsProj = verts;
                vertsProj(:,c) = 0;
%                 vertsProj
                
                
                vertBB = [min(vertsProj(:,inds(1))), max(vertsProj(:,inds(1)));
                          min(vertsProj(:,inds(2))), max(vertsProj(:,inds(2)))];
                
                u = linspace(vertBB(1,1),vertBB(1,2),nPtsPerEdgeFace);
                w = linspace(vertBB(2,1),vertBB(2,2),nPtsPerEdgeFace);
                [U,W] = meshgrid(u,w);
                in = inpolygon(U(:),W(:),vertsProj(:,inds(1)),vertsProj(:,inds(2)));

                vertsProj = [vertsProj(:,inds(1)),vertsProj(:,inds(2)); U(in), W(in)];
                vertsProj(:,3) = 0;
                vertsFinal = projectPointsToPlane(vertsProj,planeNormal);
                
                plot3(vertsFinal(:,1),vertsFinal(:,2),vertsFinal(:,3),'o')
                
%                 tri = delaunay(vertsProj(:,1),vertsProj(:,2));
                
                
%                 vertsFinal = minComp*ones(size(vertsProj,1),3);
%                 vertsFinal(:,inds(1)) =  vertsProj(:,1);
%                 vertsFinal(:,inds(2)) =  vertsProj(:,2);
%                 vertsFinal = minComp*ones(size(vertsProj,1),3);
%                 vertsFinal(:,1) =  vertsProj(:,1);
%                 vertsFinal(:,2) =  vertsProj(:,2);
%                     
%                 
%                 vertsFinal = projectPointsToPlane(vertsFinal,planeNormal);
%                                 
% %                 trisurf(tri,verts(:,1),verts(:,2),verts(:,3));
%                 trisurf(tri,vertsFinal(:,1),vertsFinal(:,2),vertsFinal(:,3));
                
                
                
                
            end
                
        case CS_CYL
            
        case CS_SPHERE
            
        otherwise
            
    end
    
end


set(gca,'NextPlot', holdStatus);