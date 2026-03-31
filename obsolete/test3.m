%% MAIN VEM 2D MIXED STOKE
close all; clearvars
addpath ./MESH/

coords = [0 0; 1 0; 0 1];
geom = get_geom(coords);

[mean_p1,mean_p2,mean_p3] = compute_mean(coords,geom);