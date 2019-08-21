%function [normals] = curvature (T,p)

%Compute normals at each vertex

nnodes = length(p);
nfaces = length(T);
weight_counter = zeros(nnodes,1);
vnormals = zeros(nnodes,3);
fnormals = zeros(nfaces,3);
coords = zeros(nnodes,3);
coords(:,:,2) = zeros(nnodes,3);
weights = zeros(nfaces,3);

T1 = T(:,1);
T2 = T(:,2);
T3 = T(:,3);



for i = 1:nnodes
	ind1 = find(T1 == i);
	ind2 = find(T2 == i);
	ind3 = find(T3 == i);
	ind = [ind1;ind2;ind3];	
	faces = T(ind,:);
	[n, e] = size(faces);
    normal = zeros(n,3);
	for j = 1:n
		edge1 = p(faces(j,2),:) - p(faces(j,1),:);
		edge2 = p(faces(j,1),:) - p(faces(j,3),:);
		normal(j,:) = (cross(edge1,edge2)/(norm(cross(edge1,edge2),2)));
    end 
	
    fnormal = sum(normal)./j; 
	unitnorm = sum(normal) ./(norm(sum(normal),2));
	vnormals(i,:) = unitnorm;
% 	
 	edgepar = dot(edge1/norm(edge1,2),vnormals(i,:))*vnormals(i,:);
% 
 	nodecoordu(i,:) = (edge1-edgepar)/norm(edge1-edgepar,2);
 	nodecoordv(i,:) = cross(vnormals(i,:), nodecoordu(i,:))/norm(cross(vnormals(i,:), nodecoordu(i,:)),2);
    randvect = rand(1,3);
	
    
	

end 

for i = 1:nfaces
	faceindex = T(i,:)';
	barycen = mean(p(faceindex,:));
	barycenter(i,:) = barycen;

	P = p(faceindex(1),:);
	Q = p(faceindex(2),:);
	R = p(faceindex(3),:);
	PQ = Q - P;
	PR = R - P;
	QR = R - Q;
	
	normf = cross(PQ,PR);
	unitnorm = normf./norm(normf,2);
	fnormals(i,:) = -unitnorm;
	facecoordu(i,:) = PQ/norm(PQ,2);
	facecoordv(i,:) = cross(fnormals(i,:), facecoordu(i,:));
	
	weights(i,1) = WeightFinder(P,Q,R,PR,PQ,QR);
	weights(i,2) = WeightFinder(Q,R,P,PQ,QR,PR);
	weights(i,3) = WeightFinder(R,Q,P,PR,QR,PQ);

end



	
faceWeits = zeros(2,2,nfaces);
vertWeits = zeros(2,2,nnodes);
faceigs = zeros(2,nfaces);
verteigs = zeros(2,nnodes);

for i = 1:nfaces
	%Here be dragons. 

	faceindex = T(i,:)';
	edge0 = p(faceindex(3),:) - p(faceindex(2),:);
	edge1 = p(faceindex(1),:) - p(faceindex(3),:);
	edge2 = p(faceindex(2),:) - p(faceindex(1),:);

	ndiff0 = vnormals(faceindex(3),:) - vnormals(faceindex(2),:);
	ndiff1 = vnormals(faceindex(1),:) - vnormals(faceindex(3),:);
	ndiff2 = vnormals(faceindex(2),:) - vnormals(faceindex(1),:);

	x1 = [dot(edge0,facecoordu(i,:));dot(edge0, facecoordv(i,:))];
	b1 = [dot(ndiff0,facecoordu(i,:));dot(ndiff0,facecoordv(i,:))];

	x2 = [dot(edge1,facecoordu(i,:));dot(edge1, facecoordv(i,:))];
	b2 = [dot(ndiff1,facecoordu(i,:));dot(ndiff1,facecoordv(i,:))];

	x3 = [dot(edge2,facecoordu(i,:));dot(edge2, facecoordv(i,:))];
	b3 = [dot(ndiff2,facecoordu(i,:));dot(ndiff2,facecoordv(i,:))];

	A = [x1' 0 ; 0 x1' ; x2' 0 ; 0 x2'; x3' 0; 0 x3'];
	b = [b1;b2;b3];


	

	linweit = A\b;


	faceWeit = [linweit(1), linweit(2); linweit(2), linweit(3)];
	
	faceWeits(:,:,i) = faceWeit;
	faceigs(:,i) = eig(faceWeit);

	

	for j = 1:3
		weight_counter(faceindex(j)) = weight_counter(faceindex(j)) + weights(i,j) ;
		ax = cross(fnormals(i,:), vnormals(faceindex(j),:));
		ax = ax / (norm(ax,2));
		theta = acos(dot(fnormals(i,:),vnormals(faceindex(j),:)));
		ux  = ax(1);
		uy = ax(2);
		uz = ax(3);
		C = cos(theta);
		S = sin(theta);
		t = 1 - C;
		rot_mat = [t*ux^2+C t*ux*uy-S*uz t*ux*uz+S*uy; t*ux*uy+S*uz t*uy^2+C t*uy*uz-S*ux; t*ux*uz-S*uy t*uy*uz+S*ux t*uz^2+C];
		u_p = rot_mat * nodecoordu(faceindex(j),:)'	;	
		v_p = rot_mat * nodecoordv(faceindex(j),:)';
		u_f = facecoordu(i,:);
		v_f = facecoordv(i,:);

		e_p = ([dot(u_p,u_f), dot(u_p,v_f)] * faceWeit * [dot(u_p,u_f); dot(u_p,v_f)]);
		f_p = ([dot(u_p,u_f), dot(u_p,v_f)] * faceWeit * [dot(v_p,u_f); dot(v_p,v_f)]);
		g_p = ([dot(v_p,u_f), dot(v_p,v_f)] * faceWeit * [dot(v_p,u_f); dot(v_p,v_f)]);
		
		vertWeits(:,:,faceindex(j)) = vertWeits(:,:,faceindex(j)) + weights(i,j)*[e_p, f_p; f_p, g_p];
	end
		

end

for i = 1:nnodes
	vertWeits(:,:,i) = vertWeits(:,:,i)/weight_counter(i);
	verteigs(:,i) = eig(vertWeits(:,:,i));
end 

	

		
