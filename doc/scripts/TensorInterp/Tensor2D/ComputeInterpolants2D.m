clc, clear all

% vectors of spatial data
x = linspace(0,1,4);
y = linspace(0,1,4);

% declare cell variables
syms u_0_0 u_0_1 u_0_2 u_1_0 u_1_1 u_1_2 u_2_0 u_2_1 u_2_2

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

% compute integrals
Q1_old = -9/256*(P(0,2) + P(1,2) + P(2,0) + P(2,1)) ...
     + 1/256*(P(0,0) + P(0,1) + P(1,0) + P(1,1) + P(2,2));

% compute u^{l+1}_{2i-1,2j-1}
Q1 = 0; Q2 = 0; Q3 = 0; Q4 = 0; Q5 = 0; Q6 = 0;
for i = 0:length(x)-1
    for j = 0:length(y)-1

        % compute P_{x_{i+1/2},y_{j+1/2}}
        Q1 = Q1 + P(i-1,j-1) * Lagrange(i,x,0.5*(x(2)+x(3))) ...
             * Lagrange(j,y,0.5*(y(2)+y(3)));

        % compute P_{x_{i},y_{j+1/2}}
        Q2 = Q2 + P(i-1,j-1) * Lagrange(i,x,x(2)) * Lagrange(j,y,0.5*(y(2)+y(3)));

        % compute P_{x_{i+1/2},y_{j}}
        Q3 = Q3 + P(i-1,j-1) * Lagrange(i,x,0.5*(x(2)+x(3))) * Lagrange(j,y,y(2));

        % compute P_{x_{i+1},y_{j+1/2}}
        Q4 = Q4 + P(i-1,j-1) * Lagrange(i,x,x(3)) * Lagrange(j,y,0.5*(y(2)+y(3)));

        % compute P_{x_{i+1/2},y_{j+1}}
        Q5 = Q5 + P(i-1,j-1) * Lagrange(i,x,0.5*(x(2)+x(3))) * Lagrange(j,y,y(3));

        % compute P_{x_{i+1},y_{j+1}}
        Q6 = Q6 + P(i-1,j-1) * Lagrange(i,x,x(3)) * Lagrange(j,y,y(3));

    end
end

% compute u^{l+1}_{2i,2j}
u4 = 4*(Q6 - Q4 - Q5 + Q1)

% compute u^{l+1}_{2i-1,2j}
u3 = 4*(Q5 - Q1 - P(0,1) + Q2)

% compute u^{l+1}_{2i-1,2j-1}
u2 = 4*(Q1 - Q2 - Q3 + P(0,0))

% compute u^{l+1}_{2i,2j-1}
u1 = 4*(Q4 - Q1 - P(1,0) + Q3)
