function [num,denom] = findFraction(decim,thresh)

    for i = 1 : 65536
        for j = 1 : 65536 
            if abs(j/i-abs(decim)) < thresh
                disp(j)
                disp(i)
            end
        end
    end

end
