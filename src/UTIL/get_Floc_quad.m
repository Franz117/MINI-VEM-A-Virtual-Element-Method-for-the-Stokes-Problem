function [F1_loc, F2_loc] = get_Floc_quad(coords, f1, f2)

    geom = get_geom(coords);
    area = geom.area;
    xv = coords(:,1);
    yv = coords(:,2);
    
    % 1. Basis functions
    x1 = coords(1,1); y1 = coords(1,2);
    x2 = coords(2,1); y2 = coords(2,2);
    x3 = coords(3,1); y3 = coords(3,2);
    x4 = coords(4,1); y4 = coords(4,2);

    J  = (x2-x1)*(y4-y1) - (y2-y1)*(x4-x1);
    s = @(x,y) ( (x-x1)*(y4-y1) - (y-y1)*(x4-x1) ) / J;
    t = @(x,y) (-(x-x1)*(y2-y1) + (y-y1)*(x2-x1) ) / J;

    p1 = @(x,y) (1 - s(x,y)).*(1 - t(x,y));
    p2 = @(x,y)      s(x,y) .*(1 - t(x,y));
    p3 = @(x,y)      s(x,y) .*      t(x,y);
    p4 = @(x,y) (1 - s(x,y)).*      t(x,y);
    %b = @(x,y) 16 *s(x,y).*(1 - s(x,y)).*t(x,y).*(1 - t(x,y));
    b = @(x,y) 8 * (1 + s(x,y) + t(x,y)) .* s(x,y).*(1-s(x,y)) .* t(x,y).*(1-t(x,y));

    
    % 2. Quadrature (Weights include area)
    Q  = get_quadrature_quad(4, coords);
    xq = Q.Points(:,1);
    yq = Q.Points(:,2);
    wq = Q.Weights;
    F1 = [f1(xq,yq)];
    F2 = [f2(xq,yq)];
    P = [p1(xq,yq),p2(xq,yq),p3(xq,yq),p4(xq,yq),b(xq,yq)];

    % 3. Build F1_loc, F2_loc (integral of f . phi_i)
    F1_loc = P' * (wq .* F1);
    F2_loc = P' * (wq .* F2);

end