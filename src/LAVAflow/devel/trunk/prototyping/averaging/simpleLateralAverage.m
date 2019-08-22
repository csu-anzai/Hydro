function [val1d, areas] = simpleLateralAverage(radFaces,phiFaces,var,doAreaNormalization)

nRadCells = numel(radFaces)-1;
nPhiCells = numel(phiFaces)-1;

val1d = zeros(nRadCells,1);
areas = zeros(nRadCells,1);
for r=1:nRadCells
    
    totalArea = 0;
    
    for p=1:nPhiCells
        
        dPhi = phiFaces(p+1)-phiFaces(p);
        
%         cellArea = dPhi/2*(radFaces(r+1)^2 - radFaces(r)^2);
        cellArea = dPhi*(radFaces(r+1)+radFaces(r))/2;
        
        val1d(r) = val1d(r) + cellArea*var(p,r);
        totalArea = totalArea + cellArea;
        
    end
    areas(r) = totalArea;
    if(doAreaNormalization)
        val1d(r) = val1d(r)/totalArea;
    end
end
        
    
