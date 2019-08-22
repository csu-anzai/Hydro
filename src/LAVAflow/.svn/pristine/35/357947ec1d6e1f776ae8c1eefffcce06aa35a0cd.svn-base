function [pi1, pia, pib] = computeProjectionIntegrals(vertsProj)

nFaceVerts = size(vertsProj,1);

pi1 = 0;
pia = 0;
pib = 0;



for i=1:nFaceVerts
    indStart = rem(i-1,nFaceVerts)+1;
    indEnd = rem(i,nFaceVerts)+1;
    vStart = vertsProj(indStart,:);
    vEnd   = vertsProj(indEnd,:);
        
    a0 = vStart(1);
    b0 = vStart(2);
    a1 = vEnd(1);
    b1 = vEnd(2);
    
    da = a1 - a0;
    db = b1 - b0;
    
    C1 = a1+a0;
    Ca = a1*C1+a0*a0;
    Cb = b1*b1 + b1*b0 + b0*b0;
    
    pi1 = pi1 + db*C1;
    pia = pia + db*Ca;
    pib = pib + da*Cb;
end

pi1 =  pi1/2;
pia =  pia/6;
pib = -pib/6;