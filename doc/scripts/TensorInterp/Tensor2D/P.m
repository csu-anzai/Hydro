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

end
