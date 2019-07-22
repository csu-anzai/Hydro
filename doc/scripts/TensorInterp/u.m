function a = u(i,j)

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

    % cases
    if (i==-1 & j==-1)
        a = u_0_0;
    elseif (i==-1 & j==0)
        a = u_0_1;
    elseif (i==-1 & j==1)
        a = u_0_2;
    elseif (i==0 & j==-1)
        a = u_1_0;
    elseif (i==0 & j==0)
        a = u_1_1;
    elseif (i==0 & j==1)
        a = u_1_2;
    elseif (i==1 & j==-1)
        a = u_2_0;
    elseif (i==1 & j==0)
        a = u_2_1;
    elseif (i==1 & j==1)
        a = u_2_2;
    end

end
