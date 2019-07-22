function rectangleNew = shrinkRectangle(IM,rectangle, nDim)


rectangleIMin = rectangle(2);
rectangleIMax = rectangle(1);
rectangleJMin = rectangle(4);
rectangleJMax = rectangle(3);
rectangleKMin = rectangle(6);
rectangleKMax = rectangle(5);

% if(nDim==2)
%     anyI = any(IM(rectangle(3):rectangle(4),rectangle(1):rectangle(2)),1);
%     anyJ = any(IM(rectangle(3):rectangle(4),rectangle(1):rectangle(2)),2);
% elseif(nDim==3)
% %     size(IM)
% %     [rectangle(3):rectangle(4)]
% %     [rectangle(1):rectangle(2)]
% %     [rectangle(5):rectangle(6)]
%     anyI = any(any(IM(rectangle(3):rectangle(4),rectangle(1):rectangle(2),rectangle(5):rectangle(6)),2),1);
%     anyJ = any(any(IM(rectangle(3):rectangle(4),rectangle(1):rectangle(2),rectangle(5):rectangle(6)),1),2);
%     anyK = any(any(IM(rectangle(3):rectangle(4),rectangle(1):rectangle(2),rectangle(5):rectangle(6)),1),3);
% end
% 
% 
% if(any(anyI))
% %     disp('In any(anyI)')
%     rectangleIMin = rectangle(1)+find(anyI,1,'first')-1;
%     rectangleIMax = rectangle(1)+find(anyI,1,'last')-1;
% end
% 
% if(any(anyJ))
% %     disp('In any(anyJ)')
%     rectangleJMin = rectangle(3)+find(anyJ,1,'first')-1;
%     rectangleJMax = rectangle(3)+find(anyJ,1,'last')-1;
% end
% if(nDim>2)
% %     disp('In any(anyK)')
%     if(any(anyK))
%         rectangleKMin = rectangle(5)+find(anyK,1,'first')-1;
%         rectangleKMax = rectangle(5)+find(anyK,1,'last')-1;
%     end
% end


for k=rectangle(5):rectangle(6)
    for j=rectangle(3):rectangle(4)
        for i=rectangle(1):rectangle(2)
            
            if(IM(j,i,k))
                
                if(k<rectangleKMin)
                    rectangleKMin = k;
                end
                if(k>rectangleKMax)
                    rectangleKMax = k;
                end
                
                
                if(j<rectangleJMin)
                    rectangleJMin = j;
                end
                if(j>rectangleJMax)
                    rectangleJMax = j;
                end
                
                
                if(i<rectangleIMin)
                    rectangleIMin = i;
                end
                if(i>rectangleIMax)
                    rectangleIMax = i;
                end
                
            end
            
        end
    end
end



rectangleNew = [rectangleIMin, rectangleIMax, rectangleJMin, rectangleJMax, rectangleKMin, rectangleKMax];