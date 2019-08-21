clear all; close all; clc;


nPoints = 40;
radius  = 1;

theta = linspace(0,2*pi,nPoints+1);
x = zeros(nPoints+2,1);
y = zeros(nPoints+2,1);

x(2:end-1) = radius*cos(theta(1:end-1));
x(1)   = x(end-1); 
x(end) = x(2);
y(2:end-1) = radius*sin(theta(1:end-1));
y(1)   = y(end-1); 
y(end) = y(2);

pos = [x,y];


curvature = zeros(nPoints,1);


for i=2:nPoints+1
    
    
    tm1 = (pos(i,:)-pos(i-1,:));
    tm1 = tm1/norm(tm1);
    
    tp1 = (pos(i+1,:)-pos(i,:));
    tp1 = tp1/norm(tp1);
    
    Lp1 = norm((pos(i+1,:)-pos(i,:))/2);
    Lm1 = norm((pos(i,:)-pos(i-1,:))/2);
    
    angle = acos(dot(tm1,tp1));
    
    curvature(i-1) = angle/(Lp1+Lm1);
    
    
end

min(curvature)
mean(curvature)
max(curvature)