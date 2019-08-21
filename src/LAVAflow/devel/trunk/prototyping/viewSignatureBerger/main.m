clear all; close all; clc;

% Read block information
blockInfoFile = '~/code/workspace/LAVAflow/build/sigberger.blkinfo';
[header,blkInfo] = hdrload(blockInfoFile);
BLK = 1;
IMIN = 2;
IMAX = 3;
JMIN = 4;
JMAX = 5;
KMIN = 6;
KMAX = 7;
PAR = 8;
CHILD1 = 9;
CHILD2 = 10;

nBlks = size(blkInfo,1);



% Read block connectivity
blockConnFile = '~/code/workspace/LAVAflow/build/sigberger.blkconn';
connList = cell(nBlks,1);
fid = fopen(blockConnFile,'r');
%discard the header line
fgetl(fid);
% Start going from here
fileLine = fgetl(fid);
while(ischar(fileLine))
    lineValues = sscanf(fileLine,'%d');
    blkID = lineValues(1)+1;
    nNeighbors = lineValues(2);
    if(nNeighbors>0)
        neighborIDs = lineValues(3:end)+1;
        connList{blkID} = neighborIDs;
    end
    
    
    fileLine = fgetl(fid);
end
    
    
fclose(fid);








%%
figure;
hold on;
axis equal

% Draw blocks
for i=1:nBlks
    
    isLeaf = blkInfo(i,CHILD1)==-1;
    if(isLeaf)
        iBnds = [blkInfo(i,IMIN),blkInfo(i,IMAX)];
        jBnds = [blkInfo(i,JMIN),blkInfo(i,JMAX)];
        kBnds = [blkInfo(i,KMIN),blkInfo(i,KMAX)];
        
        voxel([iBnds(1), jBnds(1), kBnds(1)], ...
              [iBnds(2)-iBnds(1),jBnds(2)-jBnds(1),kBnds(2)-kBnds(1)], ...
              'r', 0.8);
%         drawnow;
    end
    
end


% Draw graph connections
numLinesDrawn = 0;
for i=1:nBlks
    if(isempty(connList{i}))
        continue;
    end
    
    iBnds = [blkInfo(i,IMIN),blkInfo(i,IMAX)];
    jBnds = [blkInfo(i,JMIN),blkInfo(i,JMAX)];
    kBnds = [blkInfo(i,KMIN),blkInfo(i,KMAX)];

    blkCenter = [mean(iBnds), mean(jBnds), mean(kBnds)];
    
    neighborIDs = connList{i};
    % Loop over neighbor blocks
    for j=1:numel(neighborIDs)
        blkNeigh = neighborIDs(j);
        
        
        iBnds = [blkInfo(blkNeigh,IMIN),blkInfo(blkNeigh,IMAX)];
        jBnds = [blkInfo(blkNeigh,JMIN),blkInfo(blkNeigh,JMAX)];
        kBnds = [blkInfo(blkNeigh,KMIN),blkInfo(blkNeigh,KMAX)];
        
        
        blkCenterNeigh = [mean(iBnds), mean(jBnds), mean(kBnds)];
    
        % Draw the line
        x = [blkCenter(1); blkCenterNeigh(1)];
        y = [blkCenter(2); blkCenterNeigh(2)];
        z = [blkCenter(3); blkCenterNeigh(3)];
        plot3(x',y',z','linewidth',2);   
        
        numLinesDrawn = numLinesDrawn + 1;
        
    end
    
        drawnow;
    
end


















