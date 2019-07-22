function res = evaluateComplicatedLineIntegralSin(muStart,muEnd,etaStart,etaEnd,xiStart, xiEnd, differentialStart, differentialEnd, p,q)


res = 0;


for i=0:p
    
    
    muTerm  = (muEnd^i)*(muStart^(p-i));
    
    for j=0:q
        
        tmp = nchoosek(p,i)*nchoosek(q,j)/nchoosek(p+q,i+j);
        
        etaTerm = (etaEnd^j)*(etaStart^(q-j));
        
        intB
        
        tmp = tmp*muTerm*etaTerm;
        
        res = res + tmp;

    end
end


differential = differentialEnd-differentialStart;
res = differential*res/(p+q+1);