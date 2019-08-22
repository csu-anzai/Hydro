clear all; close all; clc;

radBnds = [1, 10];
phiBnds = [pi/4, 3*pi/4];
nPhiCells = 10;
nRadCells = 10;

[radFaces,phiFaces,RC,PC] = create_log2d_mesh(radBnds, phiBnds, nPhiCells, nRadCells);


[RF,PF] = meshgrid(radFaces,phiFaces);
XF   = RF.*cos(PF);
YF   = RF.*sin(PF);
ZF = zeros(nPhiCells+1,nRadCells+1);
XC   = RC.*cos(PC);
YC   = RC.*sin(PC);
ZC   = zeros(size(RC));

D = sin(PC)./RC;
U = RC.^2;
H = cos(PC).*exp(-RC.^2);

figure
hold on;
h = mesh(XF,YF,ZF,'edgecolor','black');
plot(XC(:),YC(:),'k.')
xlabel('x')
ylabel('y')
view(2)

% pcolor(XF,YF,ZF)
pcolor(XC,YC,H)