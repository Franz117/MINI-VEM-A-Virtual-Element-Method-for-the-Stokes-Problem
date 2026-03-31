% VIRTUAL ELEMENT VISUALIZER


% Create PDE model
clearvars; close all;
model = createpde();

% Define triangle geometry
vertices = [0 1 0; 0 0 1];  % [x; y] coordinates
edges = [1 2 3; 2 3 1];     % connectivity of triangle edges
g = decsg([2 3 vertices(1,:) vertices(2,:)]', 'T1', ('T1')');
geometryFromEdges(model, g);

% Plot geometry to verify
figure;
pdegplot(model, 'EdgeLabels', 'on');
axis equal;
title('Triangle Domain with Edge Labels');

% Generate mesh
generateMesh(model, 'Hmax', 0.05);

% Coefficients for PDE: -Δu = f, with u = g on boundary
specifyCoefficients(model, 'm', 0, 'd', 0, 'c', 1, 'a', 0, 'f', -25.9);

% Apply Dirichlet boundary condition: u = 1 - x - y
applyBoundaryCondition(model, 'dirichlet', ...
    'Edge', 1:3, ...
    'u', @(location, state) 1 - location.x - location.y);

% Solve PDE
result = solvepde(model);
u = result.NodalSolution;

% Access mesh
[p, e, t] = meshToPet(model.Mesh);

% Coordinates of triangle nodes
nodes = p';
elements = t(1:3,:)';  % node indices for each triangle

% Compute area-weighted integral of u
mean_u = 0;
total_area = 0;

for i = 1:size(elements, 1)
    % Get vertex indices and coordinates
    idx = elements(i,:);
    coords = nodes(idx, :);
    
    % Triangle area
    A = polyarea(coords(:,1), coords(:,2));
    
    % Approximate u using average at triangle's vertices
    u_avg = mean(u(idx));
    
    % Accumulate integral
    mean_u = mean_u + u_avg * A;
    total_area = total_area + A;
end

% Final mean value
mean_u = mean_u / total_area;

% Display result
fprintf('Estimated mean value of u over the triangle: %.6f\n', mean_u);


% Plot solution
figure;
pdeplot(model, 'XYData', u, 'ZData', u, 'Mesh', 'on');
title('Solution u of Poisson Equation');
xlabel('x');
ylabel('y');
