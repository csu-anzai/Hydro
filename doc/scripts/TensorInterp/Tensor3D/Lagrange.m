function w = Lagrange(j,x,xe)

    % number of points in stencil
    n = length(x);

    % offset basis from 0-index
    j = j+1;

    % compute weights
    p = 1;
    for m=1:n
        if m ~= j
            p = p * ( xe - x(m) ) / ( x(j) - x(m) );
        end
    end
    w = p;

end
