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
    a = eval(sprintf('u_%i_%i',i+1,j+1));

end
