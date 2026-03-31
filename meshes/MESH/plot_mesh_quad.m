function plot_mesh_quad(points,elements)
    hold on;
    for i = 1:size(elements, 1)
        square_nodes = points(elements(i, :), :);
        fill(square_nodes(:,1), square_nodes(:,2), 'white', 'EdgeColor', 'black');
    end
    axis equal;
    xlim([0,1]);
    ylim([0,1]);
    title('Quadrilateral Mesh with Perturbation');
    hold off;
end