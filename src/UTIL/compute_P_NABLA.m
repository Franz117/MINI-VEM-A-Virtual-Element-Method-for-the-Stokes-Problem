function [P_NABLA,G_TILDE] = compute_P_NABLA(coords,geom)

    bar = geom.bar;
    h = geom.h;
    edges = geom.edges;
    area = geom.area;
    border = geom.border;
    normals = geom.normals;

    % Basis functions
    p1 = @(x,y) 1 + 0*x + 0*y;
    p2 = @(x,y) (x - bar(1)) / h;
    p3 = @(x,y) (y - bar(2)) / h;

    % G_TILDE matrix
    G_TILDE = zeros(3, 3);
    G_TILDE(2,2) = area / (h^2);
    G_TILDE(3,3) = area / (h^2);

    % Compute P0 projection
    P0 = edges;
    P0(2:end) = P0(2:end) + edges(1:end-1);
    P0(1) = P0(1) + edges(end);
    P0 = 0.5 * P0' / border;

    % Full G matrix
    G = G_TILDE;
    G(1,:) = [P0 * p1(coords(:,1), coords(:,2)), ...
              P0 * p2(coords(:,1), coords(:,2)), ...
              P0 * p3(coords(:,1), coords(:,2))];

    % Compute B matrix
    P2 = edges .* normals(:,1);
    P2(2:end) = P2(2:end) + P2(1:end-1);
    P2(1) = P2(1) + edges(end) * normals(end,1);
    P2 = P2' / (2*h);

    P3 = edges .* normals(:,2);
    P3(2:end) = P3(2:end) + P3(1:end-1);
    P3(1) = P3(1) + edges(end) * normals(end,2);
    P3 = P3' / (2*h);

    B = [P0; P2; P3];

    % Compute P_NABLA
    P_NABLA = G \ B;
end