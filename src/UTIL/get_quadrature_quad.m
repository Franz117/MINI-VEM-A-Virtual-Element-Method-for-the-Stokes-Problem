function Q = get_quadrature_quad(degree, coords)
%GET_QUADRATURE Quadrature points and weights for a convex quadrilateral
%   degree  = desired quadrature degree
%   coords  = 4x2 array of vertices, clockwise ordering

% Manually split the quad into two triangles
% Triangle 1: [v1, v2, v3]
% Triangle 2: [v1, v3, v4]
tri_ids = [1 2 3; 1 3 4];  
v = coords;  % shorthand

Qpts = cell(2,1);
Qwts = cell(2,1);

% Get reference triangle quadrature
q_ref = quadtriangle_table(degree);

for i = 1:2
    v1 = v(tri_ids(i,1), :);
    v2 = v(tri_ids(i,2), :);
    v3 = v(tri_ids(i,3), :);

    % Jacobian matrix for triangle mapping
    J = [v2 - v1; v3 - v1]';  % 2x2

    % Transform quadrature points to physical triangle
    Qpts{i} = (J * q_ref.Points')' + v1;

    % Adjust weights
    Qwts{i} = q_ref.Weights * abs(det(J));
end

% Merge results
Q.Points  = vertcat(Qpts{:});
Q.Weights = vertcat(Qwts{:});

end