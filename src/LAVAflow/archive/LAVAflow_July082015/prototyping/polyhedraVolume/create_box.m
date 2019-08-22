function [F,E,V] = create_box(xbnds,ybnds,zbnds)

bbToVerts = [1,1,1;
            2,1,1;
            2,2,1;
            1,2,1;
            1,1,2;
            2,1,2;
            2,2,2;
            1,2,2];
        
E = [1,2;
     2,3;
     3,4;
     4,1;
     5,6;
     6,7;
     7,8;
     8,5;
     6,2;
     7,3;
     8,4;
     5,1];
          
F = [1,2,3,4;
     8,7,6,5;
     12,5,9,1;
     3,10,7,11;
     4,11,8,12;
     2,9,6,10];


bnds = [xbnds(:) ybnds(:) zbnds(:)];          
V = zeros(size(bbToVerts));
for i=1:size(V,1)
    V(i,:) = bnds(bbToVerts(i,:)+[0,2,4]);
end

              
              
              
