function plot_mesh2(mesh, filename)

    figure;
    hold on;
    axis equal;
    axis off;   % Remove axis from the plot

    % Plot elements
    for k = 1:mesh.num_elements
        nodes = mesh.connect{k};
        polygon = mesh.coords(nodes, :);
        fill(polygon(:,1), polygon(:,2), 'w', 'EdgeColor', 'k');
    end

    hold off;

    % Save the figure
    exportgraphics(gcf, filename, 'Resolution', 300);

end