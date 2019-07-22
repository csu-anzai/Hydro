clear all; 
close all; 
clc;

meshes = {};

% filename = '~/code/workspace/lavaflow/build/Tests/half_edge_mesh_0.txt';
filenameBase = '~/code/workspace/lavaflow/build/Tests/half_edge_mesh_%d.txt';

nMeshes = 0;
for i=0:2
	filename = sprintf(filenameBase, i);
	mesh = readMesh(filename);
	
	nMeshes = nMeshes+1;
	meshes{nMeshes} = mesh;
end

%%
figure
xlim([-1,2])
ylim([-1,2])
daspect([1,1,1])
hold on;

C = linspecer(nMeshes);

for i=3:3
	h = drawMesh(gca, meshes{i}, '-');
	set(h,'color',C(i,:),'linewidth',2);
end


% for f=1:mesh.nFaces
% 	edgeList = mesh.F2E{f};
% 	for e=edgeList
% 		vertList = mesh.H2V{e};
% 		plot(mesh.V(vertList,1), mesh.V(vertList,2),'k-')
% 	end
% end








