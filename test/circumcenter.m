function [x y] = circumcenter(a,b,c)

    % compute bisectors
    x1 = ( a(1) + b(1) ) / 2;
    y1 = ( a(2) + b(2) ) / 2; 
    x2 = ( a(1) + c(1) ) / 2;
    y2 = ( a(2) + c(2) ) / 2; 
      
    % compute slopes
    m1 = - ( b(1) - a(1) ) / ( b(2) - a(2) );
    m2 = - ( c(1) - a(1) ) / ( c(2) - a(2) );

    % compute circumcenter coordinates
    x = ( -m1*x1 + y1 + m2*x2 - y2 ) / ( m2 - m1 );
    y = m1*(x-x1) + y1;

end
