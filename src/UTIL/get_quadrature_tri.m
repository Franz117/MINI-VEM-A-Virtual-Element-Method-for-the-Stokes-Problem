function Q = get_quadrature_tri(degree, coords)

    % Reference quadrature rule
    q_ref = quadtriangle_table(degree);   % Points: (nq x 2), Weights: (nq x 1)

    % Triangle vertices
    v1 = coords(1,:);
    v2 = coords(2,:);
    v3 = coords(3,:);

    % Jacobian of affine map
    J = [v2 - v1; v3 - v1]';   % 2x2

    % Map quadrature points to physical triangle
    Q.Points = (J * q_ref.Points')' + v1;

    % Scale weights by area factor
    Q.Weights = q_ref.Weights * abs(det(J));

end