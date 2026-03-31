function plot_U(mesh, U1, U2,title_str)
    % Inputs:
    %   mesh - Structure with mesh.coords(:,1)=X, mesh.coords(:,2)=Y
    %   U1   - Vector field component in the x-direction (same length as mesh.coords)
    %   U2   - Vector field component in the y-direction (same length as mesh.coords)

    X = mesh.coords(:,1);
    Y = mesh.coords(:,2);

    figure;
    hold on;
    
    % Plot the vector field using quiver
    quiver(X, Y, U1, U2, 'autoscale', 'on', 'Color', 'b', 'LineWidth', 1.5);
    
    xlabel('x');
    ylabel('y');
    zlabel('Vector Components');
    title(title_str);
    grid on;
    axis equal;
    hold off;
end