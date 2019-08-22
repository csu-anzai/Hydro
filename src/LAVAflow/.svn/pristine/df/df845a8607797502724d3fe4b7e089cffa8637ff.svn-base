function [efficiency, areaMarked, areaRect] = computeEfficiency(IM, rectangle, nDim)


areaRect = (rectangle(2)-rectangle(1)+1)*(rectangle(4)-rectangle(3)+1)*(rectangle(6)-rectangle(5)+1);

if(nDim==2)
    areaMarked = sum(sum(IM(rectangle(3):rectangle(4),rectangle(1):rectangle(2))));
elseif(nDim==3)
%     size(IM)
%     [rectangle(3):rectangle(4)]
%     [rectangle(1):rectangle(2)]
%     [rectangle(5):rectangle(6)]
    areaMarked = sum(sum(sum(IM(rectangle(3):rectangle(4),rectangle(1):rectangle(2),rectangle(5):rectangle(6)))));
end

efficiency = areaMarked/areaRect;
    