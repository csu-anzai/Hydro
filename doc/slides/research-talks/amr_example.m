clc, clear all, clf

% number of points per direction
n = 128;

% grid
x = linspace(0,1,n);
[X,Y] = meshgrid(x,x);

% compute function
F = exp(-4*(X.^2+Y.^2)) .* sin(X) .* sin(Y);

% plot function
figure(1)
[C,h] = contourf(X,Y,F); hold on;
set(h,'LineColor','none');
colormap('jet');

% blocksize
bs = 16;

% overlay AMR block 1
for i = 0:bs
    plot( linspace(0,1,16), i/bs * ones(1,16), 'k', 'linewidth', 0.5 );
    plot( i/bs * ones(1,16), linspace(0,1,16) , 'k', 'linewidth', 0.5 );
end

% overlay AMR level 2 blocks
for i = 0:2*bs
    plot( linspace(0,1,16), i/(2*bs) * ones(1,16), 'k', 'linewidth', 0.5 );
    plot( i/(2*bs) * ones(1,16), linspace(0,1,16) , 'k', 'linewidth', 0.5 );
end

% overlay third AMR level blocks
for i = 0:bs
    plot( linspace(0.25,0.5,16), (0.25+i/(4*bs)) * ones(1,16), 'k', 'linewidth', 0.5 );
    plot( (0.25+i/(4*bs)) * ones(1,16), linspace(0.25,0.5,16) , 'k', 'linewidth', 0.5 );
end
for i = 0:bs
    plot( linspace(0.25,0.5,16), (0.5+i/(4*bs)) * ones(1,16), 'k', 'linewidth', 0.5 );
    plot( (0.25+i/(4*bs)) * ones(1,16), linspace(0.5,0.75,16) , 'k', 'linewidth', 0.5 );
end
