function [foundUsablePoint, splitInd, rectangleOld, rectangleNew] = findUsableSplit(IM, rectangle, possibleSplits, dir, minCells, nDim, threshold)

% Declare splitInd
splitInd = -1;

% Get the indices spanned and the number of cells in each direction
indsSpannedI = rectangle(1):rectangle(2);
indsSpannedJ = rectangle(3):rectangle(4);
indsSpannedK = rectangle(5):rectangle(6);


numCellsI = indsSpannedI(end)-indsSpannedI(1)+1;
numCellsJ = indsSpannedJ(end)-indsSpannedJ(1)+1;
numCellsK = indsSpannedK(end)-indsSpannedK(1)+1;

numCellsDir = 0;
indsSpannedDir = [];
minCellsDir = 0;
switch(dir)
    case(1)
        numCellsDir = numCellsI;
        indsSpannedDir = indsSpannedI;
        minCellsDir = minCells(1);
    case(2)
        numCellsDir = numCellsJ;
        indsSpannedDir = indsSpannedJ;
        minCellsDir = minCells(2);
    case(3)
        numCellsDir = numCellsK;
        indsSpannedDir = indsSpannedK;
        minCellsDir = minCells(3);
    otherwise
        fprintf(2,'[findUsableSplit] ERROR: Direction == %d not supported!\n',dir);
end
        
        

% Get the mean split (may not be a possible choice)
meanSplit = mean(possibleSplits);

% Get the distance from the possible points to the
% meanSplit
splitDists = abs(possibleSplits - meanSplit);

% Sort the distances in increasing order
[splitDistsSorted, orig2sort] = sort(splitDists);
sort2orig = 1:numCellsDir;
sort2orig = sort2orig(orig2sort);

% Look at all the possible split points ranked by distance
% from the mean and determine if any of them satisfy the
% constraint that Pt>=minCells and numCells-Pt+1>=minCells
% Choose the first one that satisfies this
% if none are found, just leave
rectangleOld = rectangle;
rectangleNew = rectangle;

foundUsablePoint = false;

for i=1:numel(splitDistsSorted)
    Pt = possibleSplits(sort2orig(i));

    if( (Pt >= minCellsDir) && ( numCellsDir-Pt+1 >= minCellsDir) )
        splitInd = Pt;
    else
        continue;
    end
    
    dirInd = 2*(dir-1)+1;

    % Split new rectangle in the proper direction
    rectangleNew(dirInd:dirInd+1) = [indsSpannedDir(splitInd+1), indsSpannedDir(end)];

    % Split old rectangle in the proper direction
    rectangleOld(dirInd:dirInd+1) = [indsSpannedDir(1), indsSpannedDir(splitInd)];

    % Shrink the rectangles
    rectangleNew = shrinkRectangle(IM,rectangleNew,nDim);
    
    rectangleOld = shrinkRectangle(IM,rectangleOld,nDim);

    newRectFits =   (rectangleNew(2)-rectangleNew(1))>=minCells(1) && ...
                    (rectangleNew(4)-rectangleNew(3))>=minCells(2);
            
    oldRectFits =   (rectangleOld(2)-rectangleOld(1))>=minCells(1) && ...
                    (rectangleOld(4)-rectangleOld(3))>=minCells(2);

    if(nDim==3)
        newRectFits = newRectFits && (rectangleNew(6)-rectangleNew(5))>=minCells(3);
        oldRectFits = oldRectFits && (rectangleOld(6)-rectangleOld(5))>=minCells(3);
    end
    
    [effNew, ~, ~] = computeEfficiency(IM, rectangleNew, nDim);
    [effOld, ~, ~] = computeEfficiency(IM, rectangleOld, nDim);
    
    [effNew;effOld]
    
    if(newRectFits && oldRectFits && (effNew>threshold) && (effOld>threshold))
        [effNew, effOld]
        foundUsablePoint = true;
        break;
    end              

end

if(~foundUsablePoint)
    rectangleOld = [];
    rectangleNew = [];
end
    
