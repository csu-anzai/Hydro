function drawRectangles(rectangles, nDim)

nRectangles = size(rectangles,1);

if(nDim==2)

    for i=1:nRectangles

        iBnds = [rectangles(i,1),rectangles(i,2)];
        jBnds = [rectangles(i,3),rectangles(i,4)];

        plot([iBnds(1),iBnds(2)],[jBnds(1),jBnds(1)]);
        plot([iBnds(2),iBnds(2)],[jBnds(1),jBnds(2)]);
        plot([iBnds(2),iBnds(1)],[jBnds(2),jBnds(2)]);
        plot([iBnds(1),iBnds(1)],[jBnds(2),jBnds(1)]);

    end
    
elseif(nDim==3)
    
    for i=1:nRectangles

        iBnds = [rectangles(i,1),rectangles(i,2)];
        jBnds = [rectangles(i,3),rectangles(i,4)];
        kBnds = [rectangles(i,5),rectangles(i,6)];

        voxel([iBnds(1), jBnds(1), kBnds(1)], ...
              [iBnds(2)-iBnds(1),jBnds(2)-jBnds(1),kBnds(2)-kBnds(1)], ...
              'r', 0.5);
            
        

    end
    
    
end
    
