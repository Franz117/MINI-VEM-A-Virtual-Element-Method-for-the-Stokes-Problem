%% MESHER
close all; clearvars
addpath MESH\PolyMesher\

nn = 10*2.^(10:15);  % ~10 to ~2560 elements

for i = 1:length(nn)
    n = nn(i);
%    mesh_maker("tri",  n, 0);         % unperturbed triangle mesh
%     mesh_maker("tri",  n, "auto");    % perturbed triangle mesh
     mesh_maker("quad", n, 0);         % unperturbed quad mesh
%     mesh_maker("quad", n, "auto");    % perturbed quad mesh
%     mesh_maker("poly", n);            % polygonal mesh
end
