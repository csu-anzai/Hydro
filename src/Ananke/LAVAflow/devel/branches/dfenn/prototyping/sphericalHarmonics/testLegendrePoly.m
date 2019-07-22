clear all; close all; clc;

TOL = 1e-9;

rangeL = 0:10;
xBasis = linspace(-1.0+TOL,1.0-TOL,97);


fprintf('L\tM\tx\tError\n');

for L=rangeL
    for M=0:L
        for x=xBasis
            
            res = legendrePoly(L,M,x);
            resTrue = legendre(L,x);
            resTrue = resTrue(M+1);
            
            
            error = abs((res-resTrue)/resTrue);
            
            fprintf('%i\t%i\t%f\t%E\n',L,M,x,error);
            
            if(error>1e-5)
                fprintf('Large difference found!\n');
            end
            
        end
    end
end
        
        
    