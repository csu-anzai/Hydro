function w = lagrange_derivative_weights(x,xe)

    % number of points in stencil
    n = length(x);

    % compute weights
    for j=1:n
        w(j) = 0.0;
        for i=1:n
            if i~=j
                p = 1;
                for m=1:n
                    if m==j | m==i
                    else
                        p = p * ( xe - x(m) ) / ( x(j) - x(m) );
                    end
                end
                w(j) = w(j) + p / (x(j)-x(i));
            end
        end
    end

end
