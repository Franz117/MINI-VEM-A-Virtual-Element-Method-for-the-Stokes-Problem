%obsolete TESTS

% for ie = 1:ne
%     [dof_U1, dof_U2, dof_P] = el_to_dof(ne, nv, mesh.connect{ie}, ie);
%     
%     fprintf('Element %3d:\n', ie);
%     fprintf('  DOF_U1: %s\n', mat2str(dof_U1));
%     fprintf('  DOF_U2: %s\n', mat2str(dof_U2));
%     fprintf('  DOF_P : %s\n\n', mat2str(dof_P));
% end

% all_dofs = [];
% for ie = 1:ne
%     [dof_U1, dof_U2, dof_P] = el_to_dof(ne, nv, mesh.connect{ie}, ie);
%     
%     % Concatenate all DOFs for the current element
%     all_dofs = [all_dofs; dof_U1(:); dof_U2(:); dof_P(:)];
% end

% MAKE MESH
% n=4;
% filename = "MESH/mesh_tri_"+num2str(n^2)+".txt";
% make_mesh_triang(n,filename);
% filename = "MESH/mesh_quad_"+num2str(n^2)+".txt";
% make_mesh_quad(n,filename);
% filename = "MESH/mesh_poly_"+num2str(n^2)+".txt";
% make_mesh_poly(n^2,filename);