function [A_loc, B1_loc, B2_loc, D_loc] = make_mat_loc_MINI(coords)
    % MINI element local matrix assembly: P1 + bubble for velocity, P1 for pressure

    n_dof_u = 4;  % 3 P1 + 1 bubble
    n_dof_p = 3;  % 3 P1

    area = polyarea(coords(:,1), coords(:,2));

    % Gradients of barycentric coordinates (constant over triangle)
    grads = zeros(3,2);
    for i = 1:3
        j = mod(i,3) + 1;
        k = mod(i+1,3) + 1;
        grads(i,:) = [coords(j,2) - coords(k,2), coords(k,1) - coords(j,1)] / (2*area);
    end

    A_loc = zeros(n_dof_u);
    B1_loc = zeros(n_dof_u, n_dof_p);
    B2_loc = zeros(n_dof_u, n_dof_p);

    % --- P1-P1 stiffness ---
    for i = 1:3
        for j = 1:3
            A_loc(i,j) = area * grads(i,:) * grads(j,:)';
        end
    end

    % --- Quadrature for bubble-related terms ---
    q = quadtriangle(4, 'Domain', coords);
    TR = triangulation([1 2 3], coords);
    tri_ids = ones(size(q.Points,1),1);
    lambda = cartesianToBarycentric(TR, tri_ids, q.Points);  % [Nq x 3]
    grad_lambda = grads;

    % ∇b(x) at quadrature points
    grad_bubble = 27 * ( ...
        lambda(:,2) .* lambda(:,3) * grad_lambda(1,:) + ...
        lambda(:,1) .* lambda(:,3) * grad_lambda(2,:) + ...
        lambda(:,1) .* lambda(:,2) * grad_lambda(3,:) ...
    );  % [Nq x 2]

    % --- Bubble–Bubble stiffness ---
    A_loc(4,4) = sum(q.Weights .* sum(grad_bubble.^2, 2));

    % --- Bubble–P1 coupling ---
    for i = 1:3
        grad_phi_i = grad_lambda(i,:);  % constant
        dot_vals = grad_bubble * grad_phi_i';  % [Nq x 1]
        integral = sum(q.Weights .* dot_vals);
        A_loc(i,4) = integral;
        A_loc(4,i) = integral;
    end

    % --- Divergence matrices (B1, B2) for P1 nodes ---
    for j = 1:3
        for i = 1:3
            B1_loc(i,j) = -area * grads(i,1) / 3;
            B2_loc(i,j) = -area * grads(i,2) / 3;
        end
    end

    % --- Bubble divergence contributions ---
    phi_p = lambda;  % basis functions at quadrature points
    db_dx = grad_bubble(:,1);  % ∂b/∂x
    db_dy = grad_bubble(:,2);  % ∂b/∂y

    for j = 1:3
        B1_loc(4,j) = -sum(q.Weights .* phi_p(:,j) .* db_dx);
        B2_loc(4,j) = -sum(q.Weights .* phi_p(:,j) .* db_dy);
    end

    % --- Pressure stabilization (mass-lumping) ---
    D_loc = area * ones(n_dof_p,1) / 3;

end
