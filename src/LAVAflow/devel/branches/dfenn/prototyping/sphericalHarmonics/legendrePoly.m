function Plm = legendrePoly(L, M, x)
    
    % Check that we have acceptable inputs
    if(M<0 || M>L || abs(x)>1.0)
        Plm = NaN;
        return;
    end
    
    % Compute Pmm
    Pmm = 1.0;
    
    if(M>0)
        somx2 = sqrt((1.0-x).*(1.0+x));
        fact = 1.0;
        for i=1:M
            Pmm = Pmm*(-fact*somx2);
            fact = fact + 2;
        end       
    end
    
    
    if(L==M)
        Plm = Pmm;
        return;
    else % Compute Pm(m+1)
        
        Pmmp1 = x.*(2*M+1).*Pmm;
        
        if(L == (M+1))
            Plm = Pmmp1;
            return;
        else
            for ll=M+2:L
                Pll = (x*(2*ll-1)*Pmmp1 - (ll+M-1)*Pmm)/(ll-M);
                Pmm = Pmmp1;
                Pmmp1 = Pll;
            end
            Plm = Pll;
            return;
        end
    end
    



