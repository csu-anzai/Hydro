% Compute the right-handed coordinate system permutation that 
% maximizes the z-coordinate of the normal in the new system
%
% csPerm = getCoordSystemPermutation(polyNormal)

function csPerm = getCoordSystemPermutation(polyNormal)

[~,indOrig] = max(abs(polyNormal));

indDesired = 3;

indDiff = indDesired - indOrig;

csPerm = circshift([1,2,3]',indDiff);