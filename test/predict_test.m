clc, clear all

% define a function over the domain
f = @(x,y) 0*x + 0.1*y + 1;

% reference triangle vertices and area
Tref{1,1}.vert(1,:) = [0 0];
Tref{1,1}.vert(2,:) = [1 0];
Tref{1,1}.vert(3,:) = [.5 sqrt(3)/2];
Tref{1,1}.area = sqrt(3) / 4;
Tref{1,1}.cellavg = intm2(f,Tref{1,1}.vert) / Tref{1,1}.area;

% subdivided triangles
Tref{2,1}.vert(1,:) = [0 0];
Tref{2,1}.vert(2,:) = [.5 0];
Tref{2,1}.vert(3,:) = [.25 sqrt(3)/4];
Tref{2,1}.area = sqrt(3) / 16;
Tref{2,1}.cellavg = intm2(f,Tref{2,1}.vert) / Tref{2,1}.area;

Tref{2,2}.vert(1,:) = [.5 0];
Tref{2,2}.vert(2,:) = [1 0];
Tref{2,2}.vert(3,:) = [.75 sqrt(3)/4];
Tref{2,2}.area = sqrt(3) / 16;
Tref{2,2}.cellavg = intm2(f,Tref{2,2}.vert) / Tref{2,2}.area;

Tref{2,3}.vert(1,:) = [.25 sqrt(3)/4];
Tref{2,3}.vert(2,:) = [.75 sqrt(3)/4];
Tref{2,3}.vert(3,:) = [.5 sqrt(3)/2];
Tref{2,3}.area = sqrt(3) / 16;
Tref{2,3}.cellavg = intm2(f,Tref{2,3}.vert) / Tref{2,3}.area;

Tref{2,4}.vert(1,:) = [.25 sqrt(3)/4];
Tref{2,4}.vert(2,:) = [.75 sqrt(3)/4];
Tref{2,4}.vert(3,:) = [.5 0];
Tref{2,4}.area = sqrt(3) / 16;
Tref{2,4}.cellavg = intm2(f,Tref{2,4}.vert) / Tref{2,4}.area;

% next level of subdivision of central triangle
Tref{3,1}.vert(1,:) = [.25 sqrt(3)/4];
Tref{3,1}.vert(2,:) = [.375 sqrt(3)/8];
Tref{3,1}.vert(3,:) = [.5 sqrt(3)/4];
Tref{3,1}.area = sqrt(3) / 64;
Tref{3,1}.cellavg = intm2(f,Tref{3,1}.vert) / Tref{3,1}.area;

Tref{3,2}.vert(1,:) = [.375 sqrt(3)/8];
Tref{3,2}.vert(2,:) = [.625 sqrt(3)/8];
Tref{3,2}.vert(3,:) = [.5 0];
Tref{3,2}.area = sqrt(3) / 64;
Tref{3,2}.cellavg = intm2(f,Tref{3,2}.vert) / Tref{3,2}.area;

Tref{3,3}.vert(1,:) = [.625 sqrt(3)/8];
Tref{3,3}.vert(2,:) = [.5 sqrt(3)/4];
Tref{3,3}.vert(3,:) = [.75 sqrt(3)/4];
Tref{3,3}.area = sqrt(3) / 64;
Tref{3,3}.cellavg = intm2(f,Tref{3,3}.vert) / Tref{3,3}.area;

Tref{3,4}.vert(1,:) = [.375 sqrt(3)/8];
Tref{3,4}.vert(2,:) = [.625 sqrt(3)/8];
Tref{3,4}.vert(3,:) = [.5 sqrt(3)/4];
Tref{3,4}.area = sqrt(3) / 64;
Tref{3,4}.cellavg = intm2(f,Tref{3,4}.vert) / Tref{3,4}.area;

% predict values of subdivided cell averages
Tref{3,1}.predict = Tref{2,4}.cellavg + ( Tref{2,1}.cellavg + ...
                                          Tref{2,3}.cellavg - ...
                                         2.0*Tref{2,2}.cellavg ) / 6;

Tref{3,2}.predict = Tref{2,4}.cellavg + ( Tref{2,1}.cellavg + ...
                                          Tref{2,2}.cellavg - ...
                                        2.0*Tref{2,3}.cellavg ) / 6;

Tref{3,3}.predict = Tref{2,4}.cellavg + ( Tref{2,2}.cellavg + ...
                                          Tref{2,3}.cellavg - ...
                                        2.0*Tref{2,1}.cellavg ) / 6;

Tref{3,4}.predict = Tref{2,4}.cellavg;

% compute details
Tref{3,1}.detail = abs( Tref{3,1}.predict - Tref{3,1}.cellavg );
Tref{3,2}.detail = abs( Tref{3,2}.predict - Tref{3,2}.cellavg );
Tref{3,3}.detail = abs( Tref{3,3}.predict - Tref{3,3}.cellavg );
Tref{3,4}.detail = abs( Tref{3,4}.predict - Tref{3,4}.cellavg );

Tref{3,1}.detail
Tref{3,2}.detail
Tref{3,3}.detail
Tref{3,4}.detail
