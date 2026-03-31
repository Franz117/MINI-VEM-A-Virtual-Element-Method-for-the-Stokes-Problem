function Q = get_border_quadrature_2(geom, coords, nq)
% Return quadrature points and weights along all edges of a polygon
% Inputs:
%   geom   : geometry info (bar, area, etc)
%   coords : nv x 2 coordinates of polygon vertices
%   nq     : number of Gauss points per edge
%
% Output:
%   Q.Points  : (nv*nq) x 2 array of quadrature points
%   Q.Weights : (nv*nq) x 1 array of weights along each edge
%   Q.Edges   : (nv*nq) x 1 edge index for each point
%   Q.t       : (nv*nq) x 1 local coordinate in [0,1] along edge

if nargin < 3
    nq = 2; % default 2 Gauss points per edge
end

nv = size(coords,1);
Q.Points  = zeros(nv*nq,2);
Q.Weights = zeros(nv*nq,1);
Q.Edges   = zeros(nv*nq,1);
Q.t       = zeros(nv*nq,1);

% 1D Gauss quadrature on [0,1]
[gauss_t, gauss_w] = gauss1D(nq);
% lgwt is Legendre-Gauss nodes/weights on [a,b]

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
        Q.Points(count,:) = [x1, y1] + t_local * edge_vec;
        Q.Weights(count) = gauss_w(k) * edge_len; % scale by edge length
        Q.Edges(count) = e;
        Q.t(count)     = t_local;
        count = count + 1;
    end
end

end
