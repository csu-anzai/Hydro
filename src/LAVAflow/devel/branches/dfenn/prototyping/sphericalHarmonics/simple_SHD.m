clear all; close all; clc;

L = 1;
M = 1;
nPtsTheta = 1000;
nPtsPhi   = 1000;
nPts = [nPtsTheta,nPtsPhi];
radius = 1.0;

[Ylm_1,THETA,PHI,~,~,~] = spharm(2,0,nPts,0);
[Ylm_2,~,~,~,~,~] = spharm(2,1,nPts,0);

F = Ylm_1 + Ylm_2 + spharm(2,2,nPts,0);
% F = ones(size(THETA));

X = radius*cos(THETA).*sin(PHI);
Y = radius*sin(THETA).*sin(PHI);
Z = radius*cos(PHI);

% figure
% surf(X,Y,Z,conj(Ylm))

%% Try to create E(l) = Sum_m=-l^l[|Sum[Ylm(theta,phi)*f(theta,phi)dOmega]|^2
Leval = 0:3;
E = zeros(numel(Leval),1);


for l=Leval
    mEval = -l:l;
    for m=mEval
        
        if(m>=0)
            [YlmTmp,~,~,~,~,~] = spharm(l,m,nPts,0);
        else
            [YlmTmp,~,~,~,~,~] = spharm(l,-m,nPts,0);
            YlmTmp = (-1)^abs(m)*conj(YlmTmp);
        end
        
        tmpSum = 0.0;
        
        for i=1:numel(THETA)
            
            
            tmpSum = tmpSum + conj(YlmTmp(i))*F(i);
                
        end
        
        tmpSum = abs(tmpSum^2);
        
        E(l+1) = E(l+1) + tmpSum;
        
    end
end


%%
stem(Leval,E);