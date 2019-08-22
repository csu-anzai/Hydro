function [isIntersecting, tParam, pIntersection] = linePlaneIntersection(pStart, pEnd, pPlane, normalPlane)

ERRMAX = 1e-10;

pIntersection = zeros(size(pStart));

isIntersecting = false;
tParam = NaN;
pIntersection(:) = NaN;


% make sure the points arent coincident
if(norm(pStart-pEnd)<ERRMAX)
    fprintf(2,'[Warning] Points are indistinguishable\n');
    return;
end


% Get the vector that represents the infinite line
ray = pEnd-pStart;

denom = dot(ray,normalPlane);

if(abs(denom)/norm(ray)<ERRMAX)
%     fprintf(2,'[Warning] Line lies in the plane\n');
    return;
else
    tParam = dot(pPlane-pStart, normalPlane)/denom;
    pIntersection = pStart + tParam*ray;
    if(tParam>-ERRMAX && tParam < 1+ERRMAX)
        isIntersecting = true;
    else
        isIntersecting = false;
    end
end


