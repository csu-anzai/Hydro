function res = evaluateSimpleLineIntegral(alphaStart,alphaEnd,betaStart,betaEnd,gammaStart, gammaEnd, differentialStart, differentialEnd, p,q,r)


res = 0;


for i=0:p
    
    
    alphaTerm  = (alphaEnd^i)*(alphaStart^(p-i));
    
    for j=0:q
        
        betaTerm = (betaEnd^j)*(betaStart^(q-j));
        
        for k=0:r
            
            gammaTerm = (gammaEnd^k)*(gammaStart^(r-k));
            
            tmp = nchoosek(p,i)*nchoosek(q,j)*nchoosek(r,k)/nchoosek(p+q+r,i+j+k);


            tmp = tmp*alphaTerm*betaTerm*gammaTerm;

            res = res + tmp;
            
        end

    end
end


differential = differentialEnd-differentialStart;
res = differential*res/(p+q+r+1);