function plot_point(points)
    figure;
    plot(points(:,1), points(:,2), 'bo', 'MarkerSize', 5, 'MarkerFaceColor', 'b');
    axis equal;
    xlim([0,1]);
    ylim([0,1]);
    xlabel('X');
    ylabel('Y');
    title(sprintf('Grid of Points'));
    grid on;
end