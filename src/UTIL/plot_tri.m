function plot_tri(f,xv,yv)
    % 1. Create a fine grid
    [X, Y] = meshgrid(linspace(min(xv), max(xv), 100), linspace(min(yv), max(yv), 100));
    
    % 2. Find points inside the triangle
    in = inpolygon(X, Y, xv, yv);
    
    % 3. Calculate Z and mask the outside
    Z = f(X, Y);
    Z(~in) = NaN; % Set points outside to NaN so they don't render
    
    % 4. Plot
    surf(X, Y, Z, 'EdgeColor', 'none');
