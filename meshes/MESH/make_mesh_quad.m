function mesh = make_mesh_quad(n, perturb_factor)
%MAKE_MESH_QUAD_STRUCT Generate quadrilateral mesh (optionally perturbed)
%   n: number of divisions per axis (rows = columns = n)
%   perturb_factor: max perturbation of internal nodes
%
%   Output: mesh struct with fields:
%       mesh.type
%       mesh.num_nodes
%       mesh.coords
%       mesh.num_elements
%       mesh.connect
%       mesh.boundary_nodes

%% Step 1: Create grid points
[X, Y] = meshgrid(linspace(0,1,n+1), linspace(0,1,n+1));
coords = [X(:), Y(:)];
num_nodes = size(coords,1);

%% Step 2: Create quadrilateral connectivity (CCW)
connect = {}; 
elem_id = 1;
for i = 1:n
    for j = 1:n
        node1 = (i-1)*(n+1) + j;
        node2 = node1 + 1;
        node3 = node2 + (n+1);
        node4 = node1 + (n+1);
        % store as CCW
        connect{elem_id} = [node1, node4, node3, node2];
        elem_id = elem_id + 1;
    end
end

num_elements = length(connect);

%% Step 3: Identify boundary nodes
bottom = 1:(n+1);
top    = num_nodes-(n+1)+1 : num_nodes;
left   = 1:(n+1):num_nodes;
right  = (n+1):(n+1):num_nodes;
boundary_nodes = unique([bottom, top, left, right]);

if perturb_factor > 0
    internal_nodes = setdiff(1:size(coords,1), boundary_nodes);
    coords(internal_nodes, :) = coords(internal_nodes, :) + ...
                                perturb_factor*(rand(length(internal_nodes),2));
end

%% Step 4: Create mesh struct
mesh = struct();
mesh.type = 'Quadrilateral';
mesh.num_nodes = num_nodes;
mesh.coords = coords;
mesh.num_elements = num_elements;
mesh.connect = connect;
mesh.boundary_nodes = boundary_nodes;




%% Step 5: Plot mesh
plot_mesh(mesh);
filename = sprintf('mesh_quad_p1_%d.txt', mesh.num_elements);
save_mesh(mesh, filename);

end