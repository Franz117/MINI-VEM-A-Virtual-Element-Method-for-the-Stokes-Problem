function [is_div_free, is_p_zero_mean, is_u1_zero_bc, is_u2_zero_bc] = check_prob(u1, u2, p, tol)

    if nargin < 4
        tol = 1e-6; % Default tolerance
    end

    % 1. Check divergence-free: ∂u1/∂x + ∂u2/∂y ≈ 0
    [X, Y] = meshgrid(linspace(0,1,50), linspace(0,1,50)); % grid for finite diff
    h = 1e-4;

    % Numerical derivatives
    du1_dx = (u1(X+h, Y) - u1(X-h, Y)) / (2*h);
    du2_dy = (u2(X, Y+h) - u2(X, Y-h)) / (2*h);

    div_u = du1_dx + du2_dy;
    is_div_free = max(abs(div_u(:))) < tol;

    % 2. Check mean pressure is zero (using numerical integration)
    P_vals = p(X, Y);
    area = trapz(Y(:,1), trapz(X(1,:), P_vals, 2)); % double integral
    is_p_zero_mean = abs(area) < tol;

    % 3. Check boundary conditions for u1
    y_vals = linspace(0,1,50);
    x_vals = linspace(0,1,50);

    is_u1_zero_bc = ...
        all(abs(u1(0, y_vals)) < tol) && ...
        all(abs(u1(1, y_vals)) < tol) && ...
        all(abs(u1(x_vals, 0)) < tol) && ...
        all(abs(u1(x_vals, 1)) < tol);

    % 4. Check boundary conditions for u2
    is_u2_zero_bc = ...
        all(abs(u2(0, y_vals)) < tol) && ...
        all(abs(u2(1, y_vals)) < tol) && ...
        all(abs(u2(x_vals, 0)) < tol) && ...
        all(abs(u2(x_vals, 1)) < tol);

end
