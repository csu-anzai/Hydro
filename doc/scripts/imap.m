function x = imap(l,i,j)

    iCellsInt = [4 8 16];
    jCellsInt = [1 1 1];

    x = 0;
    for m = 1:l-1
        x = x + iCellsInt(m) * jCellsInt(m);
    end

    x = x + (j-1) + jCellsInt(l) * (i-1);
    x = x+1;

end
