function [r1,c1,r2,c2,r3,c3,r4,c4,r5,c5,r6,c6] = ...
         block_dof_indices(dof_U1, dof_U2, dof_P)
%BLOCK_DOF_INDICES Generate row/column indices for block sparse assembly
% Blocks:
%  (1) U1-U1
%  (2) U2-U2
%  (3) U1-P
%  (4) U2-P
%  (5) P-U1
%  (6) P-U2

    nu = numel(dof_U1);
    np = numel(dof_P);

    % --- U1-U1 ---
    r1 = kron(dof_U1, ones(nu,1));
    c1 = kron(ones(nu,1), dof_U1);

    % --- U2-U2 ---
    r2 = kron(dof_U2, ones(nu,1));
    c2 = kron(ones(nu,1), dof_U2);

    % --- U1-P ---
    r3 = kron(dof_U1, ones(np,1));
    c3 = kron(ones(nu,1), dof_P);

    % --- U2-P ---
    r4 = kron(dof_U2, ones(np,1));
    c4 = kron(ones(nu,1), dof_P);

    % --- P-U1 ---
    r5 = kron(dof_P, ones(nu,1));
    c5 = kron(ones(np,1), dof_U1);

    % --- P-U2 ---
    r6 = kron(dof_P, ones(nu,1));
    c6 = kron(ones(np,1), dof_U2);

end