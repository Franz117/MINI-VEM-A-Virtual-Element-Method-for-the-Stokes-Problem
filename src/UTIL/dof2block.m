function [r, c] = dof2block(dof_1, dof_2)
%DOF2BLOCK Generate row/column indices for A(dof_1, dof_2)

    dof_1 = dof_1(:);
    dof_2 = dof_2(:);
    n1 = numel(dof_1);
    n2 = numel(dof_2);

    % Row indices: dof_1 cycles fastest (column-major)
    r = repmat(dof_1, n2, 1);
    % Column indices: each dof_2 repeated n1 times
    c = kron(dof_2, ones(n1,1));
end