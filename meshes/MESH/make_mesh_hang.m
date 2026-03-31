function make_mesh_hang(n)
%MAKE_MESH_HANG Generate triangular/quadrilateral mesh with hanging nodes
%   n: number of divisions per axis (rows = columns = n)
%   filename: output mesh file in the exact format you provided


%% Step 0: Parameters
x_horiz = linspace(0,1,2*n+1); % horizontal points
y_vert  = linspace(0,1,n+1);   % vertical points

%% Step 1: Create grid points
[X,Y] = meshgrid(x_horiz, y_vert);
coords = [X(:), Y(:)];
num_nodes = size(coords,1);

%% Step 2: Count elements
num_elements = (2*n+1)*n; % each row has 2*n triangles
connect = {}; 
elem_id = 1;

%% Step 1: Left side (first upward triangle of each row)
for i = 1:n
    connect{elem_id} = [i, i+1+(n+1), i+1];
    elem_id = elem_id + 1;
end
%% Step 2: Upward pyramids (corrected indices)
for j = 1:n
    while j + 2*(n+1) <= n + 2*n*(n+1)
        connect{elem_id} = [j, j+(n+1), j+2*(n+1), j+1+(n+1)];
        elem_id = elem_id + 1;
        j = j + 2*(n+1); % move to next segment (skip one)
    end
end
%% Step 3: Downward triangles (quads)
for r = 2:n+1
    j = r + (n+1);
    c=1;
    while c<= n-1
        c=c+1;
        connect{elem_id} = [j, j-1+(n+1), j+2*(n+1), j+(n+1)];
        elem_id = elem_id + 1;
        j = j + 2*(n+1); % move to next downward segment
    end
end
%% Step 4: Right-side triangles
for r = 1:n
    i = num_nodes - r; 
    connect{elem_id} = [i, i+1, i+1-(n+1)];
    elem_id = elem_id + 1;
end



%% Step 5: Identify boundary nodes
bottom = 1:(n+1);
top    = num_nodes - n : num_nodes;
left   = 1:(n+1):num_nodes;
right  = (n+1):(n+1):num_nodes;
boundary_nodes = unique([bottom, top, left, right]);

%% Step 5: Create mesh struct
mesh = struct();
mesh.type = 'Polygonal';
mesh.num_nodes = num_nodes;
mesh.coords = coords;
mesh.num_elements = num_elements;
mesh.connect = connect;
mesh.boundary_nodes = boundary_nodes;

k = 1:2:(2*n-1);          % odd numbers: 1,3,5,...,2n-1
base = k*(n+1);           % k(n+1)
remove_ids = [base(1)+1, reshape([base(2:end); base(2:end)+1],1,[])];
mesh = remove_nodes(mesh, remove_ids);
plot_mesh(mesh);
filename = sprintf('MESH/mesh_hang_%d.txt', mesh.num_elements);
save_mesh(mesh, filename);