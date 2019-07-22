function [radFaces,phiFaces,radCellCenters,phiCellCenters] = create_linear2d_mesh(radBnds, phiBnds, nPhiCells, nRadCells)

% phiMin   = pi/4;
% phiMax   = 3*pi/4;
% radMin   = 1.0;
% radMax   = 10;

phiMin = phiBnds(1);
phiMax = phiBnds(2);
radMin = radBnds(1);
radMax = radBnds(2);

radFaces = linspace(radMin,radMax,nRadCells+1);
phiFaces = linspace(phiMin,phiMax,nPhiCells+1);
% [RF,PF] = meshgrid(radFaces,phiFaces);
% XF   = RF.*cos(PF);
% YF   = RF.*sin(PF);
% ZF = zeros(nPhiCells+1,nRadCells+1);

radCellCenters = linspace(radMin + (radFaces(2)-radFaces(1))/2, radMax - (radFaces(end)-radFaces(end-1))/2,nRadCells);
phiCellCenters = linspace(phiMin + (phiFaces(2)-phiFaces(1))/2, phiMax - (phiFaces(end)-phiFaces(end-1))/2,nPhiCells);