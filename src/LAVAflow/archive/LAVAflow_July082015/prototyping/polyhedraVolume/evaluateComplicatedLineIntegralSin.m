function res = evaluateComplicatedLineIntegralSin(alphaStart,alphaEnd,betaStart,betaEnd,gammaStart, gammaEnd, xiStart, xiEnd, differentialStart, differentialEnd, p,q,r)

res = 0;

% alphaVector = [alphaStart, alphaEnd]
% betaVector = [betaStart, betaEnd]
% gammaVector = [gammaStart, gammaEnd]

for i=0:p
    
    
    alphaTerm  = (alphaEnd^i)*(alphaStart^(p-i));
    
    for j=0:q
        
        betaTerm = (betaEnd^j)*(betaStart^(q-j));
        
        for k=0:r
            
            gammaTerm = (gammaEnd^k)*(gammaStart^(r-k));
            
            tmp = nchoosek(p,i)*nchoosek(q,j)*nchoosek(r,k)/nchoosek(p+q+r,i+j+k);

            intBsinxi = integrateBSinXi(xiStart, xiEnd, p+q+r, i+j+k);

            tmp = tmp*alphaTerm*betaTerm*gammaTerm*intBsinxi;

            res = res + tmp;
            
        end

    end
end


differential = differentialEnd-differentialStart;
res = differential*res;%/(p+q+r+1);


% res = 0;
% 
% 
% for i=0:p
%     
%     
%     muTerm  = (muEnd^i)*(muStart^(p-i));
%     
%     for j=0:q
%         
%         tmp = nchoosek(p,i)*nchoosek(q,j)/nchoosek(p+q,i+j);
%         
%         etaTerm = (etaEnd^j)*(etaStart^(q-j));
%         
%         intBsinxi = integrateBSinXi(xiStart, xiEnd, p+q, i+j);
%         
%         tmp = tmp*muTerm*etaTerm*intBsinxi;
%         
%         res = res + tmp;
% 
%     end
% end
% 
% 
% differential = differentialEnd-differentialStart;
% res = differential*res/(p+q+1);



function value = integrateBSinXi(xiStart, xiEnd, n, k)
    if(n==0 && k==0)
        disp('hopefully not in here')
        h = (xiStart-xiEnd);
        if(abs(h)<1e-6)
            value = -sin(xiStart);
        else
            value = (cos(xiEnd)-cos(xiStart))/h;
        end
    else
        coeff = pi*2^(-n-3)*gamma(k+1)*gamma(n-k+1)*nchoosek(n,k);

        term1 = 4*sin(xiStart)*hypergeometricPFQRegularized([k+1,k+2]/2,[1,n+2,n+3]/2,-0.25*(xiStart-xiEnd)^2);
        term2 = (k+1)*(xiStart-xiEnd)*cos(xiStart)*hypergeometricPFQRegularized([k+2,k+3]/2,[3,n+3,n+4]/2,-0.25*(xiStart-xiEnd)^2);

        value = coeff*(term1 - term2);
    end