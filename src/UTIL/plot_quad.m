function plot_quad(f,xv,yv)
    % 1. Create a fine grid
    [X, Y] = meshgrid(linspace(min(xv), max(xv), 100), linspace(min(yv), max(yv), 100));
    
    Z = f(X, Y);
    
    % 4. Plot
    surf(X, Y, Z, 'EdgeColor', 'none');
