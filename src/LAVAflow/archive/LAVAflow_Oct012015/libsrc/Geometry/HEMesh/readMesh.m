function mesh = readMesh(filename)

fid = fopen(filename,'r');

tline = fgetl(fid);
tline = fgetl(fid);
sz = str2double(regexp(tline,'\s+','split'));

nFaces = sz(1);
nHE    = sz(2);
nVerts = sz(3);

F2E = {};
for i=1:nFaces
	tline = fgetl(fid);
	tmp = str2double(regexp(tline,'\s+','split'));
	F2E{tmp(1)} = tmp(2:end);
end
	
H2V = {};
for i=1:nHE
	tline = fgetl(fid);
	tmp = str2double(regexp(tline,'\s+','split'));
	H2V{tmp(1)} = tmp(2:end);
end

V = zeros(nVerts,3);
for i=1:nVerts
	tline = fgetl(fid);
	tmp = str2double(regexp(tline,'\s+','split'));
	V(tmp(1),:) = tmp(2:end);
end

mesh = struct('nFaces',nFaces,'nHE',nHE,'nVerts',nVerts,...
	          'F2E', {{}}, 'H2V', {{}}, 'V', [[]]);

mesh.F2E = F2E;
mesh.H2V = H2V;
mesh.V   = V;


fclose(fid);