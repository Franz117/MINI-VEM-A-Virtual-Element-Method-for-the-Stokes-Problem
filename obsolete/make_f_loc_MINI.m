function [F1_loc, F2_loc] = make_f_loc_MINI(coords, f1, f2)

    % Quadrature
    degree = 2;
    q = quadtriangle(degree, 'Domain', coords);

    % Evaluate quadrature points and weights
    X = q.Points(:,1);
    Y = q.Points(:,2);
    W = q.Weights;

    % Create a single triangle and map all points to it
    TR = triangulation([1 2 3], coords);
    tri_ids = ones(size(q.Points, 1), 1);  % replicate triangle ID for each point
    lambda = cartesianToBarycentric(TR, tri_ids, q.Points);

    % MINI basis: 3 P1 nodes + 1 bubble
    N = zeros(numel(W), 4);
    N(:,1:3) = lambda;
    N(:,4) = 27 * lambda(:,1) .* lambda(:,2) .* lambda(:,3); % bubble function

    % Evaluate RHS functions
    F1 = f1(X, Y);
    F2 = f2(X, Y);

    % Assemble local load vectors
    F1_loc = N' * (W .* F1);
    F2_loc = N' * (W .* F2);

end
