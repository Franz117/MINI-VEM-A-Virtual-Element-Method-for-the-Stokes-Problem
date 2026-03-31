clearvars;close all    
nv = 3;


Phi_mat = zeros(nv, 2*nv);
idx = 1:nv;
Phi_mat(sub2ind(size(Phi_mat), idx, 2*idx - 1)) = 1;    % phi(2*i-1) = 1
Phi_mat(sub2ind(size(Phi_mat), idx, 2*idx))     = 0.5;  % phi(2*i)   = 0.5
non_first = idx > 1;
Phi_mat(sub2ind(size(Phi_mat), idx(non_first), 2*idx(non_first)-2)) = 0.5;
Phi_mat(1, end) = 0.5;