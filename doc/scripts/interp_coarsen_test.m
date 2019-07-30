clc, clear all, clf

% function params
%f = @(x) cos(80*pi*x) .* exp(-64*x.^2);     % antiderivative
f = @(x) -tanh( (x+1/3) / (2*1e-2) ) + exp( -64^2 * (x-1/3).^2 );
I = @(x1,x2) integral(f,x1,x2);             % integral of function

% wavelet params
J = 1;                              % max level
tol = 1e-5;                        % error tolerance
jcells = @(j) 2^(j+2);              % number of points at level j

% grid params
a = -1;                              % left end
b = 1;                           % right end

% initialize matrix of cell averages and detail coefficients
for j = 0:J
    dx = (b-a) / jcells(j);
    for k = 0:jcells(j)-1
        x(j+1,k+1) = a + k*dx + dx/2;
        u(j+1,k+1) = 0;
        if j < J
            d(j+1,k+1) = 0;
        end
    end
end

% compute cell averages
for j = J:-1:0

    % width of interval
    dx = ( b - a ) / jcells(j);

    % loop through cells at level j
    for k = 0 : jcells(J)-1
        x2 = a + (k+1) * dx;
        x1 = a + k * dx;
        u(j+1,k+1) = I(x1,x2) / dx;
    end
end

% forward transform
for j = J-1:-1:0

    % compute cell width at level j
    dx = ( b - a ) / jcells(j)

    % loop through cells
    for k = 0 : jcells(j)-1

        if k == 0

            % vector of cell edges
            xvec = [ (a+0) (a+dx) (a+2*dx) (a+3*dx) ];

            % compute polynomial data
            pvec = [ 0 u(j+1,1) ( u(j+1,1) + u(j+1,2) ) ( u(j+1,1) + u(j+1,2) + u(j+1,3) ) ];

            % compute weights at cell edges for odd cell at j+1{
            w1 = lagrange_weights(xvec,xvec(1)+dx/2);
            w2 = lagrange_weights(xvec,xvec(2));

            % compute prediction of odd point at level j+1
            p1 = w1 .* pvec;
            p2 = w2 .* pvec;

        elseif k == jcells(j)-1

            % vector of cell edges
            xvec = [ (a+(k-2)*dx) (a+(k-1)*dx) (a+k*dx) (a+(k+1)*dx) ];

            % compute polynomial data
            pvec = [ 0 u(j+1,(k-2)+1) ( u(j+1,(k-2)+1) + u(j+1,(k-1)+1) )...
                        ( u(j+1,(k-2)+1) + u(j+1,(k-1)+1) + u(j+1,k+1) ) ];

            % compute weights at cell edges for odd cell at j+1{
            w1 = lagrange_weights(xvec,xvec(3)+dx/2);
            w2 = lagrange_weights(xvec,xvec(4));

            % compute prediction of odd point at level j+1
            p1 = w1 .* pvec;
            p2 = w2 .* pvec;

        else

            % vector of cell edges
            xvec = [ (a+(k-1)*dx) (a+k*dx) (a+(k+1)*dx) (a+(k+2)*dx) ];

            % compute polynomial data
            pvec = [ 0 u(j+1,(k-1)+1) ( u(j+1,(k-1)+1) + u(j+1,k+1) )...
                        ( u(j+1,(k-1)+1) + u(j+1,k+1) + u(j+1,(k+1)+1) ) ];

            % compute weights at cell edges for odd cell at j+1{
            w1 = lagrange_weights(xvec,xvec(2)+dx/2);
            w2 = lagrange_weights(xvec,xvec(3));

            % compute prediction of odd point at level j+1
            p1 = w1 .* pvec;
            p2 = w2 .* pvec;

        end

        % interpolate next grid up
        u(j+2,(2*k+1)+1) = 2 * ( sum(p2) - sum(p1) );
        u(j+2,(2*k+1)) = 2*u(j+1,k+1) - u(j+2,(2*k+1)+1);

    end
end

ubad = u;

% coarsen prediction
for k = 0 : jcells(0)-1
    ubad(1,k+1) = 0.5 * ( u(2,(2*k+1)) + u(2,(2*k+1)+1) );
end

 
% print the difference between original and cells after prediction + coarsening
abs( ubad(1,:) - u(1,:) )
