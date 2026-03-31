function [F1_loc, F2_loc] = get_Floc_tri(coords, f1, f2)

    geom = get_geom(coords);
    area = geom.area;
    xv = coords(:,1);
    yv = coords(:,2);
    
    % 1. Basis functions
    p1 = @(x,y) ((yv(2)-yv(3))*(x-xv(3)) + (xv(3)-xv(2))*(y-yv(3))) / (2*area);
    p2 = @(x,y) ((yv(3)-yv(1))*(x-xv(1)) + (xv(1)-xv(3))*(y-yv(1))) / (2*area);
    p3 = @(x,y) ((yv(1)-yv(2))*(x-xv(2)) + (xv(2)-xv(1))*(y-yv(2))) / (2*area);
    b =  @(x,y) 27 * p1(x,y) .* p2(x,y) .* p3(x,y);

    
    % 2. Quadrature (Weights include area)
    Q  = get_quadrature_tri(4, coords);
    xq = Q.Points(:,1);
    yq = Q.Points(:,2);
    wq = Q.Weights;
    F1 = [f1(xq,yq)];
    F2 = [f2(xq,yq)];
    P = [p1(xq,yq),p2(xq,yq),p3(xq,yq),b(xq,yq)];

    % 3. Build F1_loc, F2_loc (integral of f . phi_i)
    F1_loc = P' * (wq .* F1);
    F2_loc = P' * (wq .* F2);

end