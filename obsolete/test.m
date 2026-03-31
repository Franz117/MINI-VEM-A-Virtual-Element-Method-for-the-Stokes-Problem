% VIRTUAL BALL VISUALIZER

clearvars; close all;

fast_poisson_visualizer(100)


function fast_poisson_visualizer(N)
    % Fast visualizer for Laplace(u) = 1, u = 0 on boundary
    if nargin < 1
        N = 100; % Grid size (NxN)
    end

    h = 1 / (N+1);                     % Grid spacing
    f = -ones(N^2, 1);                  % Right-hand side (f = 1)
    
    % Create sparse Laplacian matrix
    e = ones(N,1);
    T = spdiags([e -4*e e], [-1 0 1], N, N);
    I = speye(N);
    L = (kron(I, T) + kron(spdiags([e e], [-1 1], N, N), I)) / h^2;

    % Solve linear system
    u = L \ f;

    % Reshape solution to 2D
    U = reshape(u, N, N);

    % Add boundary (u = 0)
    U_full = zeros(N+2, N+2);
    U_full(2:end-1, 2:end-1) = U;

    % Create meshgrid for plotting
    [X, Y] = meshgrid(0:h:1, 0:h:1);

    % Plot
    figure;
    surf(X, Y, U_full, 'EdgeColor', 'none');
    colormap jet;
    colorbar;
    title('Solution to \Delta u = -1 with u = 0 on boundary');
    xlabel('x'); ylabel('y'); zlabel('u(x,y)');
    view(3);
    axis tight;
end