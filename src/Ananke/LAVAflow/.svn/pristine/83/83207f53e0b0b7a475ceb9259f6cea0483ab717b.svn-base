function [rectangles, nRectangles, efficiencies] = computeRectangleClusters(IM, minCells, threshold)

nDim = ndims(IM);

nI = 1;
nJ = 1;
nK = 1;

if(nDim==2)
    nI = size(IM,2);
    nJ = size(IM,1);
elseif(nDim==3)
    nI = size(IM,2);
    nJ = size(IM,1);
    nK = size(IM,3);
end

nRectangles = 1;
rectangles = [1,nI,1,nJ, 1, nK];

while(true)
% for iter=1:10
%     fprintf('New iteration\n');
    didSomething = false;
    
    for r=1:nRectangles
        
        % Shrink rectangles to be bounded by flagged values
        rectangles(r,:) = shrinkRectangle(IM,rectangles(r,:),nDim);
        
        % Compute block efficiency
        [efficiency, ~, ~] = computeEfficiency(IM, rectangles(r,:), nDim);
        efficiency
        if(efficiency<threshold)
            
            % Compute signatures
            [sigmaI, sigmaJ, sigmaK, diffI, diffJ, diffK] = computeSignatures(IM,rectangles(r,:), nDim);
            
            signChangedI = (diffI(1:end-1)>0 & diffI(2:end)<0) | (diffI(1:end-1)<0 & diffI(2:end)>0);
            signChangedJ = (diffJ(1:end-1)>0 & diffJ(2:end)<0) | (diffJ(1:end-1)<0 & diffJ(2:end)>0);
            
            signChangedK = [];
            if(nDim==3)
                signChangedK = (diffK(1:end-1)>0 & diffK(2:end)<0) | (diffK(1:end-1)<0 & diffK(2:end)>0);
            end
            
            foundUsablePoint = false;
                    
            rectangleOld = [];
            rectangleNew = [];
                
            % Check I-zeros
            zerosI = find(sigmaI==0);
            if(~isempty(zerosI) && ~foundUsablePoint)
                possibleSplits = zerosI;
                [foundUsablePoint, splitInd, rectangleOld, rectangleNew] = findUsableSplit(IM,rectangles(r,:), possibleSplits, 1, minCells, nDim, threshold);
            end
            
            % Check J-zeros
            zerosJ = find(sigmaJ==0);
            if(~isempty(zerosJ) && ~foundUsablePoint)
                disp('Splitting on J zeros')
                possibleSplits = zerosJ;
                [foundUsablePoint, splitInd, rectangleOld, rectangleNew] = findUsableSplit(IM,rectangles(r,:), possibleSplits, 2, minCells, nDim, threshold);
            end
            
            % Check K-zeros
            if(nDim==3)
                zerosK = find(sigmaK==0);
                if(~isempty(zerosK) && ~foundUsablePoint)
                    possibleSplits = zerosK;
                    [foundUsablePoint, splitInd, rectangleOld, rectangleNew] = findUsableSplit(IM,rectangles(r,:), possibleSplits, 3, minCells, nDim, threshold);
                end
            end
            
            % Check if there are any inflection points in I
            if(any(signChangedI) && ~foundUsablePoint)
                possibleSplits = find(signChangedI)+1;
                [foundUsablePoint, splitInd, rectangleOld, rectangleNew] = findUsableSplit(IM,rectangles(r,:), possibleSplits, 1, minCells, nDim, threshold);
            end
            
            % Check if there are any inflection points in J
            if(any(signChangedJ) && ~foundUsablePoint)
                possibleSplits = find(signChangedJ)+1;
                [foundUsablePoint, splitInd, rectangleOld, rectangleNew] = findUsableSplit(IM,rectangles(r,:), possibleSplits, 2, minCells, nDim, threshold);
            end
            
            % Check if there are any inflection points in K
            if(nDim==3)
                if(any(signChangedK) && ~foundUsablePoint)
                    possibleSplits = find(signChangedK)+1;
                    [foundUsablePoint, splitInd, rectangleOld, rectangleNew] = findUsableSplit(IM,rectangles(r,:), possibleSplits, 3, minCells, nDim, threshold);
                end
            end
                        
            
            % If we found a usable point, split the rectangle up 
            if(foundUsablePoint)
                rectangles(r,:) = rectangleOld;
                rectangles(nRectangles+1,:) = rectangleNew;
                nRectangles = nRectangles + 1;
                didSomething = true;
                break;
            end

                       
        end
    end
    
    
    
    if(~didSomething)
        break;
    end
    
    
end

efficiencies = zeros(nRectangles,1);

for r=1:nRectangles
    % Compute block efficiency
    [efficiency, ~, ~] = computeEfficiency(IM, rectangles(r,:), nDim);
    efficiencies(r) = efficiency;
end

