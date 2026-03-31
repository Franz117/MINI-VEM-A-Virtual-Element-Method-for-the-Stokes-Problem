function Q = get_border_quadrature_simple(coords, nq)
% Returns a quadrature along the polygon boundary
% Inputs:
%   coords : nv x 2 polygon vertex coordinates
%   nq     : number of Gauss points per edge (default 2)
% Outputs:
%   Q.Points  : (nv*nq) x 2 array of quadrature points
%   Q.Weights : (nv*nq) x 1 array of weights (scalar, edge length * Gauss weight)

if nargin < 2
    nq = 2;  % default 2 points per edge
end

nv = size(coords,1);
Q.Points  = zeros(nv*nq,2);
Q.Weights = zeros(nv*nq,1);

% Hard-coded 1D Gauss points on [0,1] for nq = 1,2,3
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
    % Edge nodes
    n1 = e;
    n2 = mod(e,nv)+1;

    a = coords(n1,:);
    b = coords(n2,:);
    edge_vec = b - a;
    edge_len = norm(edge_vec);

    for k = 1:nq
        t = gauss_t(k);                % mapped to [0,1]
        Q.Points(count,:)  = (1-t)*a + t*b;
        Q.Weights(count)   = gauss_w(k) * edge_len;
        count = count + 1;
    end
end

end
