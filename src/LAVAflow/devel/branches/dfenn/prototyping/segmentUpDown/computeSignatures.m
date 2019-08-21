function [sigmaI, sigmaJ, sigmaK, diffI, diffJ, diffK] = computeSignatures(IM, rectangle, nDim)

rangeI = rectangle(1:2);
rangeJ = rectangle(3:4);
rangeK = rectangle(5:6);

nI = rangeI(2)-rangeI(1)+1;
nJ = rangeJ(2)-rangeJ(1)+1;
nK = rangeK(2)-rangeK(1)+1;

sigmaI = zeros(nI,1);
sigmaJ = zeros(nJ,1);
diffI = zeros(nI-2,1);
diffJ = zeros(nJ-2,1);

sigmaK = [];
diffK = [];
if(nDim==3)
    sigmaK = zeros(nK,1);
    diffK = zeros(nK-2,1);
end
    

if(nDim==2)
    if(nI>3)
        for i=1:nI
            sigmaI(i) = sum(IM(rangeJ(1):rangeJ(2),i+rangeI(1)-1));
        end

        for i=2:nI-1
            diffI(i-1) = sigmaI(i+1) - 2*sigmaI(i) + sigmaI(i-1);
        end
    end

    if(nJ>3)
        for j=1:nJ
            sigmaJ(j) = sum(IM(j+rangeJ(1)-1,rangeI(1):rangeI(2)));
        end

        for j=2:nJ-1
            diffJ(j-1) = sigmaJ(j+1) - 2*sigmaJ(j) + sigmaJ(j-1);
        end
    end
elseif(nDim==3)
    if(nI>3)
        for i=1:nI
            sigmaI(i) = sum(sum(IM(rangeJ(1):rangeJ(2),i+rangeI(1)-1,rangeK(1):rangeK(2))));
        end

        for i=2:nI-1
            diffI(i-1) = sigmaI(i+1) - 2*sigmaI(i) + sigmaI(i-1);
        end
    end

    if(nJ>3)
        for j=1:nJ
            sigmaJ(j) = sum(sum(IM(j+rangeJ(1)-1,rangeI(1):rangeI(2),rangeK(1):rangeK(2))));
        end

        for j=2:nJ-1
            diffJ(j-1) = sigmaJ(j+1) - 2*sigmaJ(j) + sigmaJ(j-1);
        end
    end
    
    
    if(nK>3)
        for k=1:nK
            sigmaK(k) = sum(sum(IM(rangeJ(1):rangeJ(2),rangeI(1):rangeI(2),k+rangeK(1)-1)));
        end

        for k=2:nK-1
            diffK(k-1) = sigmaK(k+1) - 2*sigmaK(k) + sigmaK(k-1);
        end
    end
end
    
    
    
    
    