function Q = get_border_quadrature_AI(geom, coords)
    % 2-point Gauss quadrature on each edge
    n_edges = size(geom.edges, 1);
    n_nodes = size(coords, 1);
    
    % Gauss points (in [-1,1] reference interval)
    gp = [-1/sqrt(3), 1/sqrt(3)];
    weights = [1, 1];  % weights for reference interval [-1,1]
    
    % Initialize
    Q.Points = zeros(2*n_edges, 2);
    Q.Weights1 = zeros(2*n_edges, 1);
    Q.Weights2 = zeros(2*n_edges, 1);
    
    for e = 1:n_edges
        n1 = e;
        n2 = mod(e, n_nodes) + 1;
        
        % Edge endpoints
        a = coords(n1,:);
        b = coords(n2,:);
        
        % Map Gauss points to physical edge
        for q = 1:2
            idx = 2*(e-1)+q;
            t = 0.5*(1 + gp(q));  % map to [0,1]
            Q.Points(idx,:) = (1-t)*a + t*b;
            
            % Jacobian factor (edge length/2)
            jac = geom.edges(e)/2;
            
            Q.Weights1(idx) = weights(q) * jac * geom.normals(e,1);
            Q.Weights2(idx) = weights(q) * jac * geom.normals(e,2);
        end
    end
end