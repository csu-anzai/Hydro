function w = lagrange_weights(x,xe)

    % number of points in stencil
    n = length(x);

    % compute weights
    for j=1:n
        p = 1;
        for m=1:n
            if m ~= j
                p = p * ( xe - x(m) ) / ( x(j) - x(m) );
            end
        end
        w(j) = p;
    end

end
