clear all; close all; clc


r = linspace(0,1,100);
offset1 = 0.5;
offset2 = 0.75;

rNew1 =  100*(r-offset1);
rNew2 =  100*(r-offset2);

y = (tanh(rNew1)+1)/2 - (tanh(rNew2)+1)/2;

figure
plot(r,y);
