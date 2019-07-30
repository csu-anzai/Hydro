function a = u(i,j,k)

    % global vars
    global u_0_0_0
    global u_0_1_0
    global u_0_2_0
    global u_1_0_0
    global u_1_1_0
    global u_1_2_0
    global u_2_0_0
    global u_2_1_0
    global u_2_2_0
    global u_0_0_1
    global u_0_1_1
    global u_0_2_1
    global u_1_0_1
    global u_1_1_1
    global u_1_2_1
    global u_2_0_1
    global u_2_1_1
    global u_2_2_1
    global u_0_0_2
    global u_0_1_2
    global u_0_2_2
    global u_1_0_2
    global u_1_1_2
    global u_1_2_2
    global u_2_0_2
    global u_2_1_2
    global u_2_2_2

    % cases
    a = eval(sprintf('u_%i_%i_%i',i+1,j+1,k+1));

end
