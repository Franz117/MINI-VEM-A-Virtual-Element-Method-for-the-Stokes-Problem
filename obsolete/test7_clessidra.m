%% Define vertices and boundary conditions
vertices = [0 1 1 0; 0 0 1 1]; % [x; y]
vertex_values = [1 -1 1 -1];    % alternating values

% Triangulate the polygon
dt = delaunayTriangulation(vertices');
triangles = dt.ConnectivityList;

% Create PDE model
model = createpde();

% Define geometry
geometryFromMesh(model, vertices, triangles');

% Generate mesh
generateMesh(model, 'Hmax', 0.05);

% Laplace equation
specifyCoefficients(model, 'm', 0, 'd', 0, 'c', 1, 'a', 0, 'f', 0);

% Get mesh nodes
nodes = model.Mesh.Nodes';

% Identify boundary nodes
TR = triangulation(triangles, vertices');
boundaryEdges = freeBoundary(TR);
boundaryNodes = unique(boundaryEdges(:));

% Assign boundary values to nodes by interpolating along polygon edges
u_bc = zeros(size(nodes, 1), 1);
numVertices = size(vertices, 2);
for k = 1:length(boundaryNodes)
    n = boundaryNodes(k);
    px = nodes(n, 1); py = nodes(n, 2);
    % Find closest edge of polygon
    minDist = Inf;
    for i = 1:numVertices
        v1 = i; v2 = mod(i, numVertices) + 1;
        x1 = vertices(1, v1); y1 = vertices(2, v1);
        x2 = vertices(1, v2); y2 = vertices(2, v2);
        % Project point onto edge
        t = ((px - x1) * (x2 - x1) + (py - y1) * (y2 - y1)) / ((x2 - x1)^2 + (y2 - y1)^2);
        t = max(0, min(1, t));
        projx = x1 + t * (x2 - x1);
        projy = y1 + t * (y2 - y1);
        dist = sqrt((px - projx)^2 + (py - projy)^2);
        if dist < minDist
            minDist = dist;
            u_bc(n) = (1 - t) * vertex_values(v1) + t * vertex_values(v2);
        end
    end
end

% Apply Dirichlet BC on boundary nodes
applyBoundaryCondition(model, 'dirichlet', 'Edge', 1:model.Geometry.NumEdges, 'u', @(location, ~) interpolateBC(location, vertices, vertex_values));

% Solve PDE
result = solvepde(model);
u = result.NodalSolution;

% Plot solution
figure;
pdeplot(model, 'XYData', u, 'ZData', u, 'Mesh', 'on');
title('Harmonic function on quadrilateral');
xlabel('x');
ylabel('y');

%% Compute area-weighted mean of u
nodes = model.Mesh.Nodes';
elements = model.Mesh.Elements';
mean_u = 0;
total_area = 0;

for i = 1:size(elements, 1)
    idx = elements(i, :);
    coords = nodes(idx, :);
    tri_area = polyarea(coords(:, 1), coords(:, 2));
    u_avg = mean(u(idx));
    mean_u = mean_u + u_avg * tri_area;
    total_area = total_area + tri_area;
end

mean_u = mean_u / total_area;
fprintf('Estimated mean: %.6f\n', mean_u);

%% Helper function to interpolate boundary conditions
function val = interpolateBC(location, vertices, vertex_values)
    x = location.x;
    y = location.y;
    val = zeros(size(x));
    numVertices = size(vertices, 2);
    for k = 1:numel(x)
        px = x(k);
        py = y(k);
        minDist = Inf;
        for i = 1:numVertices
            v1 = i;
            v2 = mod(i, numVertices) + 1;
            x1 = vertices(1, v1);
            y1 = vertices(2, v1);
            x2 = vertices(1, v2);
            y2 = vertices(2, v2);
            t = ((px - x1) * (x2 - x1) + (py - y1) * (y2 - y1)) / ((x2 - x1)^2 + (y2 - y1)^2);
            t = max(0, min(1, t));
            projx = x1 + t * (x2 - x1);
            projy = y1 + t * (y2 - y1);
            dist = sqrt((px - projx)^2 + (py - projy)^2);
            if dist < minDist
                minDist = dist;
                val(k) = (1 - t) * vertex_values(v1) + t * vertex_values(v2);
            end
        end
    end
end
