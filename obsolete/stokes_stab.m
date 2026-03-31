%% MAIN VEM 2D MIXED STOKE
close all; clearvars
addpath(genpath('./MESH'))
addpath(genpath('./UTIL'))

% read mesh
%meshfile = "mesh_tri_p_32.txt";
%meshfile = "mesh_hang_465.txt";
%meshfile = "mesh_dia_1600.txt";
%meshfile = "mesh_quad_2601.txt";
meshfile = "mesh_poly_1280.txt";
%meshfile = "mesh_tri.txt";
mesh = read_mesh(meshfile);
%plot_mesh(mesh);

% read prob
sol = get_prob();
u1 = sol.u1;
u2 = sol.u2;
p =  sol.p;
f1 = sol.f1;
f2 = sol.f2;

% dimensionalities
ne = mesh.num_elements;
nv = mesh.num_nodes;
dim_U = nv + ne;
dim_P = nv;
ndofs = 2*dim_U + dim_P +1;

% memory alloc
row_cell = cell(ne, 1);
col_cell = cell(ne, 1);
val_cell = cell(ne, 1);
RHS_row_cell = cell(ne, 1);
RHS_val_cell = cell(ne, 1);
hh = zeros(ne,1);

% main loop
parfor ie=1:ne 
    
    % main logic
    [A_loc,B1_loc,B2_loc, D_loc,h] = get_MATloc_stab(mesh.coords(mesh.connect{ie},:));
    [F1_loc,F2_loc] = get_Floc_vem(mesh.coords(mesh.connect{ie},:),f1,f2);

    % memory management
    rows = [];
    cols = [];
    vals = [];
    [dof_U1,dof_U2,dof_P] = el_to_dof(ne,nv,mesh.connect{ie},ie);
    blocks = {
        dof_U1, dof_U1, A_loc
        dof_U2, dof_U2, A_loc
        dof_U1, dof_P , B1_loc
        dof_U2, dof_P , B2_loc
        dof_P , dof_U1, B1_loc.'
        dof_P , dof_U2, B2_loc.'
        ndofs , dof_P, D_loc
        dof_P , ndofs, D_loc.'
    };
    for k = 1:size(blocks,1)
        dof_1 = blocks{k,1}(:);
        dof_2 = blocks{k,2}(:);
        values  = blocks{k,3};
        [r, c] = dof2block(dof_1, dof_2);
        rows = [rows; r];
        cols = [cols; c];
        vals = [vals; values(:)];
    end
    row_cell{ie} = rows;
    col_cell{ie} = cols;
    val_cell{ie} = vals;
    RHS_row_cell{ie} = [dof_U1(:); dof_U2(:)];
    RHS_val_cell{ie} = [F1_loc(:);F2_loc(:)];

    % get h
    hh(ie) = h;

end
h = max(hh);

% assemble A matrix and RHS
A_row = vertcat(row_cell{:});
A_col = vertcat(col_cell{:});
A_val = vertcat(val_cell{:});
A = sparse(A_row, A_col, A_val);

F_row = vertcat(RHS_row_cell{:});
F_val = vertcat(RHS_val_cell{:});
F = accumarray(F_row, F_val, [ndofs, 1]);

% border condition (Hiron's trick)
bd_nodes = mesh.boundary_nodes;
coords_bd = mesh.coords(bd_nodes,:);

u1_bd = u1(coords_bd(:,1), coords_bd(:,2));
u2_bd = u2(coords_bd(:,1), coords_bd(:,2));

[dof_U1_bd,dof_U2_bd,~] = node_to_dof(bd_nodes,ne,nv);
all_bd = [dof_U1_bd;dof_U2_bd];

A(all_bd,:)=0;
A(sub2ind(size(A), all_bd, all_bd)) = 1;
F(dof_U1_bd) = u1_bd;
F(dof_U2_bd) = u2_bd;

% system solver
SOL = A\F;
U1 = SOL(1:dim_U);
U2 = SOL(dim_U+1:2*dim_U);
P = SOL(2*dim_U+1:end-1);

% plot
X = mesh.coords(:,1);
Y = mesh.coords(:,2);
plot_P(mesh,P,"Discrete P",true);
%plot_P(mesh,p(X,Y),"True P",true);
%plot_U(mesh, U1(1:nv), U2(1:nv),"Discrete U");
%plot_U(mesh,u1(X,Y),u2(X,Y),"True U");


% error estimation
err0_u= zeros(ne, 1);
err1_u = zeros(ne, 1);
err0_p = zeros(ne, 1);
parfor ie=1:ne

    verts = mesh.connect{ie};
    coords = mesh.coords(verts,:);
    [err0_u_loc, err1_u_loc, err0_p_loc] = get_ERRloc_vem(coords,verts,U1,U2,P,sol);

    err0_u(ie) = err0_u_loc;
    err1_u(ie) = err1_u_loc;
    err0_p(ie) = err0_p_loc;

end
err0_u = sqrt(sum(err0_u));
err1_u = sqrt(sum(err1_u));
err0_p = sqrt(sum(err0_p));

% terminal prints
fprintf('h = %12.4e\n',h);
fprintf('Err L2 of p = %12.4e\n',err0_p);
fprintf('Err L2 of u = %12.4e\n',err0_u);
fprintf('Err H1 of u = %12.4e\n',err1_u);
