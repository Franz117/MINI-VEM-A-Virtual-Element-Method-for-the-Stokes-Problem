%% TEST
close all; clearvars
addpath ./MESH/
addpath ./quadtriangle/

meshfile = "mesh_tri_256.txt";
mesh = read_mesh(meshfile);
[u1,u2,p,f1,f2,mu] = get_prob();

% ASSEMBLY MATRIX
ne = mesh.num_elements;
nv = mesh.num_nodes;
dim_U = nv + ne;
dim_P = nv;
A = zeros(2*dim_U + dim_P +1);

for ie=1:ne

    [dof_U1,dof_U2,dof_P] = el_to_dof(ne,nv,mesh.connect{ie},ie);

    [B1_loc,B2_loc] = make_border_loc(mesh.coords(mesh.connect{ie},:));
    A(dof_U1,dof_P) = A(dof_U1,dof_P) + B1_loc;
    A(dof_U2,dof_P) = A(dof_U2,dof_P) + B2_loc;
    A(dof_P,dof_U1) = A(dof_P,dof_U1) + B1_loc';
    A(dof_P,dof_U2) = A(dof_P,dof_U2) + B2_loc';

end

%APPLY BORDER CONDITION (trucco IRON)
for ib=1:numel(mesh.boundary_nodes)
    node = mesh.boundary_nodes(ib);
    [dof_U1,dof_U2,dof_P] = node_to_dof(node,ne,nv);
    A(dof_U1,:) = 0; A(dof_U1,dof_U1) = 1;
    A(dof_U2,:) = 0; A(dof_U2,dof_U2) = 1;
    B(dof_U1)= 0; B(dof_U2)= 0;
end

spy(A)
