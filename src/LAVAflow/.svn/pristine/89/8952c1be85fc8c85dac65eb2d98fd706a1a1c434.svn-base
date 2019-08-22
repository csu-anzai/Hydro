function sum = hypergeometric2F3(a,b,z)

NMAX = 100;
ERRMAX = 1e-6;

sum = 1;

aOld = ones(size(a));
bOld = ones(size(b));
zOld = 1;
nFactOld = 1;
for n=1:NMAX
    
    aNew = (a+n-1).*aOld;
    bNew = (b+n-1).*bOld;
    zNew = zOld*z;
    nFactNew = nFactOld*n;
    
    tmp = prod(aNew)/prod(bNew)*(zNew/nFactNew);
    
    sumOld = sum;
    sum = sum + tmp;
    error = abs((sum-sumOld)/sum);
    if(error<ERRMAX)
        fprintf('2F3 Converged!\n\tn = %d\n\terror = %5.5E\n',n,error);
        return;
    end
    
    aOld = aNew;
    bOld = bNew;
    zOld = zNew;
    nFactOld = nFactNew;
end
        

fprintf(2,'Failed to converge!\n\tn = %d\n\terror = %5.5E\n',n,error);
