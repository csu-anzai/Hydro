function [F1,Fa,Fb,Fg] = computeFaceIntegrals(vertsTrans, polyNormalTrans)



w = -dot(polyNormalTrans,vertsTrans(1,:));

% vertsProj = projectPointsToPlane(vertsTrans,polyNormalTrans);
% vertsProj = projectPointsToPlane(vertsTrans,[0,0,1]);
vertsProj = vertsTrans;
vertsProj(:,3) = 0;

% w = -dot(polyNormalTrans,vertsProj(1,:));

k1 = 1/polyNormalTrans(3);
k2 = k1*k1;

[pi1, pia, pib] = computeProjectionIntegrals(vertsProj);

F1 = pi1;

Fa = k1*pia;
Fb = k1*pib;
Fg = -k2*(polyNormalTrans(1)*pia + polyNormalTrans(2)*pib + w*pi1);