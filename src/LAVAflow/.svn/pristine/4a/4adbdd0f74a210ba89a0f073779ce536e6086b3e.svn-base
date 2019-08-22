function res = evaluateComplicatedLineIntegralCos(muStart,muEnd,etaStart,etaEnd,xiStart, xiEnd, differentialStart, differentialEnd, p,q)


res = 0;


for i=0:p
    
    
    muTerm  = (muEnd^i)*(muStart^(p-i));
    
    for j=0:q
        
        tmp = nchoosek(p,i)*nchoosek(q,j)/nchoosek(p+q,i+j);
        
        etaTerm = (etaEnd^j)*(etaStart^(q-j));
        
        intBcosxi = integrateBCosXi(xiStart, xiEnd, p+q, i+j);
        
        tmp = tmp*muTerm*etaTerm*intBcosxi;
        
        res = res + tmp;

    end
end


differential = differentialEnd-differentialStart;
res = differential*res/(p+q+1);



function value = integrateBCosXi(xiStart, xiEnd, n, k)

    if(n==0 && k==0)
        h = (xiStart-xiEnd);
        if(abs(h)<1e-6)
            value = cos(xiStart);
        else
            value = (sin(xiStart)-sin(xiEnd))/h;
        end
    else
        coeff = pi*2^(-n-3)*gamma(k+1)*gamma(n-k+1)*nchoosek(n,k);

        term1 = 4*cos(xiStart)*hypergeometricPFQRegularized([k+1,k+2]/2,[1,n+2,n+3]/2,-0.25*(xiStart-xiEnd)^2);
        term2 = (k+1)*(xiStart-xiEnd)*sin(xiStart)*hypergeometricPFQRegularized([k+2,k+3]/2,[3,n+3,n+4]/2,-0.25*(xiStart-xiEnd)^2);

        value = coeff*(term1 + term2);
    end