clear all; close all; clc;

radBnds = [1, 2];
phiBnds = [0,pi/3];
nPhiCells = 100;
nRadCells = 10;

% [radCellFaces,phiCellFaces,radCellCenters,phiCellCenters] = create_log2d_mesh(radBnds, phiBnds, nPhiCells, nRadCells);
[radCellFaces,phiCellFaces,radCellCenters,phiCellCenters] = create_linear2d_mesh(radBnds, phiBnds, nPhiCells, nRadCells);
[RC,PC] = meshgrid(radCellCenters,phiCellCenters);

D = RC;
U = sin(PC);
H = cos(PC);
DbarExact = @(r) r;
UbarExact = @(r) (cos(phiBnds(1))-cos(phiBnds(2)))/(phiBnds(2)-phiBnds(1))*ones(size(r));
HbarExact = @(r) (sin(phiBnds(2))-sin(phiBnds(1)))/(phiBnds(2)-phiBnds(1))*ones(size(r));
HprimeExact = @(R,T) cos(T) - HbarExact(R);

FcKernelExact = @(R,T) R.*sin(T).*(cos(T) - HbarExact(R));


Dbar = simpleLateralAverage(radCellFaces,phiCellFaces,D,true);
Ubar = simpleLateralAverage(radCellFaces,phiCellFaces,U,true);
Hbar = simpleLateralAverage(radCellFaces,phiCellFaces,H,true);

fprintf('Norm of Hbar: %f\n',norm(Hbar-HbarExact(radCellCenters')));

Hprime = zeros(size(H));
for r=1:nRadCells
    Hprime(:,r) = H(:,r) - Hbar(r);
end

% Hprime = HprimeExact(RC,PC);

FcKernel = D.*U.*Hprime;
% FcKernel = FcKernelExact(RC,PC);

[Fc,areas] = simpleLateralAverage(radCellFaces,phiCellFaces,FcKernel,false);
Fc = Fc;



%%
% DbarExact = @(r) (-cos(phiBnds(2))+cos(phiBnds(1)))/(phiBnds(2)-phiBnds(1))./r;
% UbarExact = @(r) r.^2;
% HbarExact = @(r) (sin(phiBnds(2))-sin(phiBnds(1)))/(phiBnds(2)-phiBnds(1))*exp(-r.^2);
% HprimeExact = @(R,T) exp(-R.^2).*cos(T) - HbarExact(R);

norm(Hprime - HprimeExact(RC,PC))

DbarErr = (Dbar-DbarExact(radCellCenters'))./DbarExact(radCellCenters');
UbarErr = (Ubar-UbarExact(radCellCenters'))./UbarExact(radCellCenters');
HbarErr = (Hbar-HbarExact(radCellCenters'))./HbarExact(radCellCenters');

figure
% subplot(1,2,1);
plot(radCellCenters,DbarErr, ...
     radCellCenters,UbarErr, ...
     radCellCenters,HbarErr);
legend('density','velocity','enthalpy')
    
%% sdfsdf
% FcCoeff = -(cos(phiBnds(2))-cos(phiBnds(1)))*((phiBnds(2)-phiBnds(1))*(cos(phiBnds(2))+cos(phiBnds(1)))-2*sin(phiBnds(2))+2*sin(phiBnds(1)))/(2*(phiBnds(2)-phiBnds(1)));
% FcExact = @(r) FcCoeff*r.^2.*exp(-r.^2);
% FcErr = (Fc-FcExact(radCellCenters'))./FcExact(radCellCenters');

FcCoeff = -cos(2*phiBnds(2)) + cos(2*phiBnds(1)) + 2*(sin(2*phiBnds(2))+sin(2*phiBnds(1)) - 2*sin(phiBnds(2)-phiBnds(1)))/(phiBnds(2)-phiBnds(1));
FcExact = @(r) 0.25*FcCoeff*r.^2;
FcErr = (Fc-FcExact(radCellCenters'))./FcExact(radCellCenters');
% figure
% plot(radCellCenters,FcErr);

%%
figure
hold on;
plot(radCellCenters,FcExact(radCellCenters'))
plot(radCellCenters,Fc,'k-o')

%%
[RF,PF] = meshgrid(radCellFaces,phiCellFaces);
XF   = RF.*cos(PF);
YF   = RF.*sin(PF);
ZF = zeros(nPhiCells+1,nRadCells+1);
XC   = RC.*cos(PC);
YC   = RC.*sin(PC);
ZC   = zeros(size(RC));

figure
hold on;
h = mesh(XF,YF,ZF,'edgecolor','black','edgealpha',0);
plot(XC(:),YC(:),'k.')
xlabel('x')
ylabel('y')
view(2)

% pcolor(XF,YF,ZF)
% hand=pcolor(XC,YC,H);
hand=pcolor(XC,YC,Hprime);
set(hand,'edgealpha',0)