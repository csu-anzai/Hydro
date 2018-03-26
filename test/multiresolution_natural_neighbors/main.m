clc, clear all, clf

% number of grid points in each direction
numPts = 20;

% max level
J = 4;
ptMin = numPts^2 / J;

% detail threshold
threshold = 1e-2;

% bounding box
a = -4;
b = 4;

% define gaussian functional on domain
f = @(x,y) 1 * exp( -( .5*x.^2 + .5*y.^2) );

% create linear arrays
x = linspace(a,b,numPts);
y = linspace(a,b,numPts);

% create grid array
[X,Y] = meshgrid(x,y);

% add perturbation to interior points
eps = .75e-1;
for i=2:length(X)-1
    for j = 2:length(Y)-1
        X(i,j) = X(i,j) + 2*eps*rand - eps;
        Y(i,j) = Y(i,j) + 2*eps*rand - eps;
    end
end

% create suitable array for matlab built-in function
k = 1;
for i=1:length(X)
    for j=1:length(X)
        P(k,1) = X(i,j);
        P(k,2) = Y(i,j);
        k = k + 1;
    end
end

% decompose mesh
for j = J:-1:1
    
    % compute the delaunay triangulation
    DT{j} = delaunayTriangulation(P);

    % compute number of odd nodes
    numOdd = floor( length(P) / 2 );

    % select odd vertices randomly
    oddNodes = 1 + floor( length(P)*rand(1,numOdd) );
    
    % ensure all odd vertices have enough neighbors
    for i = 1:numOdd
        neighbors = neighborFind( DT{j}, oddNodes(i) );
        if length( neighbors ) > 2
        else
            flag = false;
            while flag == false
                newNode = 1 + floor( length(P)*rand );
                if ismember(oddNodes,newNode) == zeros(1,length(oddNodes)) ...
                        & length( neighborFind( DT{j}, newNode ) ) > 2
                    oddNodes(i) = newNode;
                    flag = true;
                end
            end
        end
    end

    % compute the detail coefficients
    for i = 1:numOdd

        % recompute 1-ring neighbors
        neighbors = neighborFind( DT{j}, oddNodes(i) );
     
        % compute the node's approximate functional value via natural neighbors
        vq = griddata( P(neighbors,1), P(neighbors,2), f(P(neighbors,1), P(neighbors,2)),...
                P(oddNodes(i),1), P(oddNodes(i),2) );

        % compute actual value based on known function
        vtrue = f( P(oddNodes(i),1), P(oddNodes(i),2) );

        % compute the detail coefficient
        detailVec(i) = vtrue - vq;

    end

    % delete nodes below threshold 
    k = 1;
    for i = 1:numOdd
        if abs( detailVec(i) ) < threshold
            deleteNodes(k) = oddNodes(i);
            k = k + 1;
        end
    end
    if k > 1
        P(deleteNodes,:) = [ ];
        clear deleteNodes;
    end
    clear oddNodes;
      
end

% plot the triangulation
figure(1)
triplot(DT{J});

figure(2)
triplot(DT{1});

% plot the function
figure(3)
surf(X,Y,f(X,Y));
