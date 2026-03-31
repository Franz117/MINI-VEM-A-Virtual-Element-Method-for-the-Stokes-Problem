function Q = get_quadrature(degree,coords) %faster

    % for non convex
    edges = [(1:size(coords,1))', [2:size(coords,1) 1]'];
    DT = delaunayTriangulation(coords, edges);
    inside = isInterior(DT);

    tri_ids = DT.ConnectivityList(inside,:);                % (nTriangles x 3)
    points = DT.Points;                                     % (nPoints x 2)
    all_tricoords = reshape(points(tri_ids',:), 3, [], 2);  % (3 x nTriangles x 2)
    all_tricoords = permute(all_tricoords, [2 1 3]);        % (nTriangles x 3 x 2)

    nt = size(all_tricoords,1);
    Qpts = cell(nt, 1);
    Qwts = cell(nt, 1);

    q_ref = quadtriangle_table(degree);
    for i=1:nt
        v1 = [all_tricoords(i,1,1) all_tricoords(i,1,2) ];
        v2 = [all_tricoords(i,2,1) all_tricoords(i,2,2) ];
        v3 = [all_tricoords(i,3,1) all_tricoords(i,3,2) ];
        J = [v2 - v1; v3 - v1]';

        Qpts{i} = (J * q_ref.Points')' + v1;   % Transform each point
        Qwts{i} = q_ref.Weights * abs(det(J)); % Adjust weights with absolute value of determinant of Jacobian
    end
    Q.Points = vertcat(Qpts{:});
    Q.Weights = vertcat(Qwts{:});

return

%% OLD
% function Q = get_quadrature(degree, coords)
% 
%     % Construct edges of the polygon
%     edges = [(1:size(coords,1))', [2:size(coords,1), 1]'];
% 
%     % Use constrained Delaunay triangulation
%     DT = delaunayTriangulation(coords, edges);
% 
%     % Keep only triangles strictly inside the polygon
%     inside = isInterior(DT); % <-- This filters out triangles outside the non-convex region
%     T = DT.ConnectivityList(inside, :); % <-- Filtered triangle connectivity
% 
%     nt = size(T, 1);
%     Qpts = cell(nt, 1);
%     Qwts = cell(nt, 1);
%     for i = 1:nt
%         tricoords = DT.Points(T(i,:), :);
%         q = quadtriangle(degree, 'Domain', tricoords);
%         Qpts{i} = q.Points;
%         Qwts{i} = q.Weights;
%     end
%     Q.Points = vertcat(Qpts{:});
%     Q.Weights = vertcat(Qwts{:});
% 
% return