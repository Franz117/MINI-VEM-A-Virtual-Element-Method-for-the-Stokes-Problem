function plot_mesh(mesh)

    figure;
    hold on;
    axis equal;
    title(sprintf('%s Mesh', mesh.type));  % Dynamic title
    xlabel('X');
    ylabel('Y');

    % Plot elements
    for k = 1:mesh.num_elements
        nodes = mesh.connect{k};
        polygon = mesh.coords(nodes, :);
        fill(polygon(:,1), polygon(:,2), 'w', 'EdgeColor', 'k');
    end

    % Highlight boundary nodes (optional)
    if isfield(mesh, 'boundary_nodes')
        plot(mesh.coords(mesh.boundary_nodes,1), ...
             mesh.coords(mesh.boundary_nodes,2), ...
             'ro', 'MarkerFaceColor', 'r');
    end

    % Optional: number the nodes
    for i = 1:mesh.num_nodes
        text(mesh.coords(i,1), mesh.coords(i,2), ...
             sprintf('%d', i), ...
             'Color', 'b', 'FontSize', 8, 'HorizontalAlignment', 'center');
    end

    hold off;
end