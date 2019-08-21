clear all; close all; clc;

% Create image
nI = 101;
nJ = 101;
diameter = nI^2+nJ^2;

IM = zeros(nI,nJ);

% U shape
rInner = diameter/8^2;
rOuter = diameter/4^2;
[I,J] = meshgrid(1:nI,1:nJ); R = (I-floor(nI/2)).^2+(J-floor(nJ/2)).^2;
IM = R>=rInner & R<=rOuter;% & J<=floor(nJ/2);


% Circle near the top-right
origin = [51,51];
R = (I-origin(1)).^2+(J-origin(2)).^2;
IM = IM | R<=100;


la = IM;

figure
hold on;
% subplot(2,2,1);
plot(I(la(:)),J(la(:)),'x')
xlim([0,nI+1])
ylim([0,nJ+1])
axis equal

%% Compute segmentation


threshold = 0.8;

minCellsI = 3;
minCellsJ = 3;
minCellsK = 1;
minCells = [minCellsI, minCellsJ, minCellsK];

[rectangles, nRectangles] = computeRectangleClusters(IM, minCells, threshold);


nRectangles 

drawRectangles(rectangles,2);







