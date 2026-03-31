function perturbed_points = perturb(points, epsilon)
    % Identify the boundary points
    x = points(:,1);
    y = points(:,2);
    
    isBoundary = (x == 0 | x == 1 | y == 0 | y == 1);
    
    perturbation = (rand(size(points)) - 0.5) * 2 * epsilon;

    perturbed_points = points;
    perturbed_points(~isBoundary, :) = perturbed_points(~isBoundary, :) + perturbation(~isBoundary, :);
    
    perturbed_points = max(0, min(1, perturbed_points));
end