clear all; close all; clc;

xbnds = [0.5,1];
ybnds = [1,2];
zbnds = [-3,3];

[F,E,V] = create_box(xbnds,ybnds,zbnds);

figure;

draw_box(gca,F,E,V,'cartesian')

%%
% xbnds = [0.5,1];
% ybnds = [pi/4,3*pi/4];
% zbnds = [-1,2];
% 
% [F,E,V] = create_box(xbnds,ybnds,zbnds);
% 
% figure;
% 
% draw_box(gca,F,E,V,'cylindrical')

%%
% xbnds = [0.5,1];
% ybnds = [pi/4,3*pi/4];
% zbnds = [0,pi/2];
% 
% [F,E,V] = create_box(xbnds,ybnds,zbnds);
% 
% figure;
% 
% draw_box(gca,F,E,V,'spherical')