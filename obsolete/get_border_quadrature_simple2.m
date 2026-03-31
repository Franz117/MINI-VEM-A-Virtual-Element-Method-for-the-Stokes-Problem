function Q = get_border_quadrature_simple2(coords, nq)
% Return a simple quadrature over the polygon boundary
% Inputs:
%   coords : nv x 2 coordinates of polygon vertices
%   nq     : number of Gauss points per edge (default 2)
%
% Output:
%   Q.Points  : (nv*nq) x 2 array of quadrature points
%   Q.Weights : (nv*nq) x 1 array of weights

if nargin < 2
    nq = 2; % default 2 Gauss points per edge
end

nv = size(coords,1);
Q.Points  = zeros(nv*nq,2);
Q.Weights = zeros(nv*nq,1);

% Simple hard-coded 1D Gauss points on [0,1]
switch nq
    case 1
        gauss_t = 0.5;
        gauss_w = 1;
    case 2
        gauss_t = [0.211324865405187, 0.788675134594813]';
        gauss_w = [0.5, 0.5]';
    case 3
        gauss_t = [0.112701665379258, 0.5, 0.887298334620742]';
        gauss_w = [5/18, 8/18, 5/18]';
    otherwise
        error('nq>3 not implemented');
end

count = 1;
for e = 1:nv
    % edge nodes
    n1 = e;
    n2 = mod(e,nv)+1;

    x1 = coords(n1,1); y1 = coords(n1,2);
    x2 = coords(n2,1); y2 = coords(n2,2);

    edge_vec = [x2-x1, y2-y1];
    edge_len = norm(edge_vec);

    for k = 1:nq
        t_local = gauss_t(k);
        Q.Points(count,:)  = [x1, y1] + t_local*edge_vec;
        Q.Weights(count)   = gauss_w(k) * edge_len;  % scale by edge length
        count = count + 1;
    end
end

end