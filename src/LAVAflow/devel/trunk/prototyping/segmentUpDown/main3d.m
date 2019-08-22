clear all; close all; clc;

% Create image
nI = 100;
nJ = 100;
nK = 100;
diameter = sqrt(nI^2+nJ^2 + nK^2);

IM = zeros(nI,nJ,nK);

% U shape
rInner = diameter/5;
rOuter = diameter/4;
[I,J,K] = meshgrid(1:nI,1:nJ,1:nK); R = sqrt((I-floor(nI/2)).^2+(J-floor(nJ/2)).^2+(K-floor(nK/2)).^2);
IM = R>=rInner & R<=rOuter & J<=floor(nJ/2);


% % Circle near the top-right
origin = [80,80,80];
R = sqrt((I-origin(1)).^2+(J-origin(2)).^2 + (K-origin(3)).^2);
IM = IM | R<=15;

la = IM;
nSkip = 4;
laNew = la(1:nSkip:end, 1:nSkip:end, 1:nSkip:end);

figure
hold on;
% subplot(2,2,1);
plot3(I(la),J(la),K(la),'x','markersize',5)
xlim([0,nI+1])
ylim([0,nJ+1])
zlim([0,nK+1])
axis equal

%% Compute segmentation

nDim = ndims(IM);

threshold = 0.9;

minCellsI = 5;
minCellsJ = 5;
minCellsK = 5;
minCells = [minCellsI, minCellsJ, minCellsK];

[rectangles, nRectangles, efficiencies] = computeRectangleClusters(IM, minCells, threshold);

nRectangles 

drawRectangles(rectangles, nDim);
xlabel('x')
ylabel('y')
zlabel('z')






