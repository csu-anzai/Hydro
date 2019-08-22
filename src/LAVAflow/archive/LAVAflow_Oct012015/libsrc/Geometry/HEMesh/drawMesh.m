function h = drawMesh(ax, mesh, linespec)

holdOld = get(ax,'NextPlot');
set(ax,'NextPlot','add');

h = [];

for f=1:mesh.nFaces
	if(mod(f,2)==1)
		linestyle = '-';
	else
		linestyle = '--';
	end
	
	edgeList = mesh.F2E{f};
	for e=edgeList
		vertList = mesh.H2V{e};
		h(end+1) = plot(ax, mesh.V(vertList,1), mesh.V(vertList,2),linestyle);
	end
end

set(ax,'NextPlot',holdOld);