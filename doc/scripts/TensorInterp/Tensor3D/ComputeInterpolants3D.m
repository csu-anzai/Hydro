clc, clear all

% vectors of spatial data
x = linspace(0,1,4);
y = linspace(0,1,4);
z = linspace(0,1,4);

% declare cell variables
syms u_0_0_0 u_0_1_0 u_0_2_0 u_1_0_0 u_1_1_0 u_1_2_0 u_2_0_0 u_2_1_0 u_2_2_0
syms u_0_0_1 u_0_1_1 u_0_2_1 u_1_0_1 u_1_1_1 u_1_2_1 u_2_0_1 u_2_1_1 u_2_2_1
syms u_0_0_2 u_0_1_2 u_0_2_2 u_1_0_2 u_1_1_2 u_1_2_2 u_2_0_2 u_2_1_2 u_2_2_2

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

% compute some needed integral values
Q0 = 0; Q1 = 0; Q2 = 0; Q3 = 0; Q4 = 0; Q5 = 0; Q6 = 0; Q7 = 0;
for i = 0:length(x)-1
    for j = 0:length(y)-1
        for k = 0:length(z)-1

            % compute P_{x_{i+1/2},y_{j+1/2},z_{k+1/2}}
            Q0 = Q0 + P(i-1,j-1,k-1) * Lagrange(i,x,0.5*(x(2)+x(3))) ...
                 * Lagrange(j,y,0.5*(y(2)+y(3))) ...
                 * Lagrange(k,z,0.5*(z(2)+z(3)));

            % compute P_{x_{i+1/2},y_{j+1/2},z_{k}}
            Q1 = Q1 + P(i-1,j-1,k-1) * Lagrange(i,x,0.5*(x(2)+x(3))) ...
                 * Lagrange(j,y,0.5*(y(2)+y(3))) ...
                 * Lagrange(k,z,z(2));

            % compute P_{x_{i},y_{j+1/2},z_{k+1/2}}
            Q2 = Q2 + P(i-1,j-1,k-1) * Lagrange(i,x,x(2)) ...
                 * Lagrange(j,y,0.5*(y(2)+y(3))) ...
                 * Lagrange(k,z,0.5*(z(2)+z(3)));

            % compute P_{x_{i+1/2},y_{j},z_{k+1/2}}
            Q3 = Q3 + P(i-1,j-1,k-1) * Lagrange(i,x,0.5*(x(2)+x(3))) ...
                 * Lagrange(j,y,y(2)) ...
                 * Lagrange(k,z,0.5*(z(2)+z(3)));

            % compute P_{x_{i},y_{j+1/2},z_{k}}
            Q4 = Q4 + P(i-1,j-1,k-1) * Lagrange(i,x,x(2)) ...
                 * Lagrange(j,y,0.5*(y(2)+y(3))) ...
                 * Lagrange(k,z,z(2));

            % compute P_{x_{i+1/2},y_{j},z_{k}}
            Q5 = Q5 + P(i-1,j-1,k-1) * Lagrange(i,x,0.5*(x(2)+x(3))) ...
                 * Lagrange(j,y,y(2)) ...
                 * Lagrange(k,z,z(2));

            % compute P_{x_{i},y_{j},z_{k+1/2}}
            Q6 = Q6 + P(i-1,j-1,k-1) * Lagrange(i,x,x(2)) ...
                 * Lagrange(j,y,y(2)) ...
                 * Lagrange(k,z,0.5*(z(2)+z(3)));

        end
    end
end

% compute u_{2i-1,2j-1,2k-1}
u1 = 8*(Q0 - Q1 - Q2 - Q3 + Q4 + Q5 + Q6 - P(0,0,0))
