% function plot_P(mesh, P, title_str, use_surface)
%     % Inputs:
%     %   mesh        - Structure with mesh.coords(:,1)=X, mesh.coords(:,2)=Y
%     %   P           - Vector of function values at each (X,Y)
%     %   use_surface - (Optional) true to plot surface, false to plot scatter (default: false)
% 
%     if nargin < 4
%         use_surface = false;
%     end
% 
%     X = mesh.coords(:,1);
%     Y = mesh.coords(:,2);
% 
%     figure;
%     hold on;
% 
%     if use_surface
%         n_elem = mesh.num_elements;
%         patch_handles = gobjects(n_elem,1);
% 
%         for ie = 1: n_elem
%             dof = mesh.connect{ie};
%             x = X(dof);
%             y = Y(dof);
%             p = P(dof);
%             % c = mean(u); % if solid color, else vertex coloring
%             patch_handles(ie) = patch(x, y, p, p, 'EdgeColor', 'k');
%         end
%     else
%         scatter3(X, Y, P, 50, P, 'filled');
%     end
% 
%     xlabel('x');
%     ylabel('y');
%     zlabel('p(x, y)');
%     title(title_str);
%     colorbar;
%     view(3);
%     grid on;
%     %axis equal;
%     hold off;
% end


%% OLD
function plot_P(mesh, P,title_str,use_surface)
    % Inputs:
    %   mesh        - Structure with mesh.coords(:,1)=X, mesh.coords(:,2)=Y
    %   P           - Vector of function values at each (X,Y)
    %   use_surface - (Optional) true to plot surface, false to plot scatter (default: false)

    if nargin < 4
        use_surface = false;
    end

    X = mesh.coords(:,1);
    Y = mesh.coords(:,2);

    figure;
    hold on;

    if use_surface
        % Create a grid and interpolate values to get a surface
        xlin = linspace(min(X), max(X), 100);
        ylin = linspace(min(Y), max(Y), 100);
        [Xq, Yq] = meshgrid(xlin, ylin);
        Pq = griddata(X, Y, P, Xq, Yq, 'natural'); % 'natural' or 'cubic'

        surf(Xq, Yq, Pq, 'EdgeColor', 'none'); 
        shading interp;
    else
        scatter3(X, Y, P, 50, P, 'filled');
    end

    xlabel('x');
    ylabel('y');
    zlabel('p(x, y)');
    title(title_str);
    colorbar;
    view(3);
    grid on;
    hold off;
end
