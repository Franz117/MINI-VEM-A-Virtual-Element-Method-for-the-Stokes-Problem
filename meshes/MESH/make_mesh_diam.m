function make_mesh_diam(Nx,Ny)

%% -------------------------------------------------
% Domain
%% -------------------------------------------------
Lx = 1;
Ly = 2;

hx = Lx / Nx;
hy = Ly / Ny;

%% -------------------------------------------------
% 1) Corner nodes (ROW-WISE numbering)
%% -------------------------------------------------
num_corner = (Nx+1)*(Ny+1);
coords_corner = zeros(num_corner,2);

k = 1;
for j = 0:Ny
    for i = 0:Nx
        coords_corner(k,:) = [i*hx , j*hy];
        k = k + 1;
    end
end

%% -------------------------------------------------
% 2) Center nodes (one per quad)
%% -------------------------------------------------
num_centers = Nx*Ny;
coords_center = zeros(num_centers,2);

k = 1;
for j = 0:Ny-1
    for i = 0:Nx-1
        coords_center(k,:) = [(i+0.5)*hx , (j+0.5)*hy];
        k = k + 1;
    end
end

%% Combine all coordinates
coords = [coords_corner ; coords_center];
num_nodes = size(coords,1);

%% -------------------------------------------------
% 3) Connectivity (4 triangles per quad)
%% -------------------------------------------------
connect = cell(4*Nx*Ny,1);
elem = 1;

for j = 1:Ny
    for i = 1:Nx
        
        % Corner nodes (row-wise numbering)
        n1 = (j-1)*(Nx+1) + i;      % bottom-left
        n2 = n1 + 1;                % bottom-right
        n4 = n1 + (Nx+1);           % top-left
        n3 = n4 + 1;                % top-right
        
        % Center node
        nc = num_corner + (j-1)*Nx + i;
        
        % 4 triangles (counter-clockwise)
        connect{elem} = [n1 n2 nc]; elem = elem + 1;
        connect{elem} = [n2 n3 nc]; elem = elem + 1;
        connect{elem} = [n3 n4 nc]; elem = elem + 1;
        connect{elem} = [n4 n1 nc]; elem = elem + 1;
        
    end
end

num_elements = length(connect);

%% -------------------------------------------------
% 4) Boundary nodes
%% -------------------------------------------------
tol = 1e-12;

boundary_nodes = find( ...
       abs(coords(:,1)) < tol | ...
       abs(coords(:,1)-Lx) < tol | ...
       abs(coords(:,2)) < tol | ...
       abs(coords(:,2)-Ly) < tol );

%% -------------------------------------------------
% 5) Build mesh structure
%% -------------------------------------------------
mesh = struct();

mesh.type = 'Triangular';
mesh.num_nodes = num_nodes;
mesh.coords = coords;
mesh.coords(:,2) = mesh.coords(:,2)/2;
mesh.num_elements = num_elements;
mesh.connect = connect;
mesh.boundary_nodes = boundary_nodes;

disp(mesh)
%plot_mesh(mesh)


filename = sprintf('MESH/diam/mesh_diam_%d.txt', mesh.num_elements);
save_mesh(mesh, filename);