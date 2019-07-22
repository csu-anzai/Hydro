function q = P(i,j)

    % global vars
    global u_0_0
    global u_0_1
    global u_0_2
    global u_1_0
    global u_1_1
    global u_1_2
    global u_2_0
    global u_2_1
    global u_2_2

    % return value
    q = 0;
    for ii = -1:i-1
        for jj = -1:j-1
            q = q + u(ii,jj);
        end
    end

    if (i==-1 & j==-1)
        a = 0;
    elseif (i==0 & j==-1)
        a = 0;
    elseif (i==1 & j==-1)
        a = 0;
    elseif (i==2 & j==-1)
        a = 0;
    elseif (i==-1 & j==0)
        a = 0;
    elseif (i==0 & j==0)
        a = u(-1,-1);
    elseif (i==1 & j==0)
        a = u(-1,-1) + u(0,-1);
    elseif (i==2 & j==0)
        a = u(-1,-1) + u(0,-1) + u(1,-1);
    elseif (i==-1 & j==1)
        a = 0;
    elseif (i==0 & j==1)
        a = u(-1,-1) + u(-1,0);
    elseif (i==1 & j==1)
        a = u(-1,-1) + u(-1,0) + u(0,-1) + u(0,0);
    elseif (i==2 & j==1)
        a = u(-1,-1) + u(-1,0) + u(0,-1) + u(0,0) + u(1,-1) + u(1,0);
    elseif (i==-1 & j==2)
        a = 0;
    elseif (i==0 & j==2)
        a = u(-1,-1) + u(-1,0) + u(-1,1);
    elseif (i==1 & j==2)
        a = u(-1,-1) + u(-1,0) + u(-1,1) + u(0,-1) + u(0,0) + u(0,1);
    elseif (i==2 & j==2)
        a = u(-1,-1) + u(-1,0) + u(-1,1) + u(0,-1) + u(0,0) + u(0,1) ...
            + u(1,-1) + u(1,0) + u(1,1);
    end

    disp(q-a)

end
