function vec = neighborFind(DT,index)

    k = 1;
    vec = zeros(3,1);
    for i = 1:length(DT.ConnectivityList)
        p = ismember( DT.ConnectivityList(i,:) , index );
        if p(1) == 1 | p(2) == 1 | p(3) == 1
            for j = 1:length(p)
                if p(j) == 0 & ismember( vec, DT.ConnectivityList(i,j) ) == 0
                    vec(k) = DT.ConnectivityList(i,j);
                    k = k + 1;
                end
            end
        end
    end

end
