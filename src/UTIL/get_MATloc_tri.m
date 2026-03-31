function [K_loc,B1_loc,B2_loc,D_loc,h] = get_MATloc_tri(coords)

    geom = get_geom(coords);
    h = geom.h;
    area = geom.area;
    xv = coords(:,1);
    yv = coords(:,2);
    
    % 1. Basis functions
    p1 = @(x,y) ((yv(2)-yv(3))*(x-xv(3)) + (xv(3)-xv(2))*(y-yv(3))) / (2*area);
    p2 = @(x,y) ((yv(3)-yv(1))*(x-xv(1)) + (xv(1)-xv(3))*(y-yv(1))) / (2*area);
    p3 = @(x,y) ((yv(1)-yv(2))*(x-xv(2)) + (xv(2)-xv(1))*(y-yv(2))) / (2*area);
    b =  @(x,y) 27 * p1(x,y) .* p2(x,y) .* p3(x,y);

    p1x = @(x,y) 0*x + 0*y + (yv(2)-yv(3))/ (2*area);
    p1y = @(x,y) 0*x + 0*y + (xv(3)-xv(2))/ (2*area);
    p2x = @(x,y) 0*x + 0*y + (yv(3)-yv(1))/ (2*area);
    p2y = @(x,y) 0*x + 0*y + (xv(1)-xv(3))/ (2*area);
    p3x = @(x,y) 0*x + 0*y + (yv(1)-yv(2))/ (2*area);
    p3y = @(x,y) 0*x + 0*y + (xv(2)-xv(1))/ (2*area);
    bx = @(x,y)27 * ((p2(x,y) .* p3(x,y)).* p1x(x,y) + ...
                    (p1(x,y) .* p3(x,y)).* p2x(x,y) + ...
                    (p1(x,y) .* p2(x,y)).* p3x(x,y) );
    by = @(x,y)27 * ((p2(x,y) .* p3(x,y)).* p1y(x,y) + ...
                    (p1(x,y) .* p3(x,y)).* p2y(x,y) + ...
                    (p1(x,y) .* p2(x,y)).* p3y(x,y) );
    
    % 3. Quadrature (Weights include area)
    Q  = get_quadrature_tri(4, coords);
    xq = Q.Points(:,1);
    yq = Q.Points(:,2);
    wq = Q.Weights;
    P = [p1(xq,yq),p2(xq,yq),p3(xq,yq)];
    PX = [p1x(xq,yq),p2x(xq,yq),p3x(xq,yq),bx(xq,yq)];
    PY = [p1y(xq,yq),p2y(xq,yq),p3y(xq,yq),by(xq,yq)];


    % 4. Build K_loc (Stiffness matrix: integral of grad_phi_i . grad_phi_j)
    K_loc = PX' * (wq .* PX) + PY' * (wq .* PY);

    % 5. Build B1_loc, B2_loc (Divergence: integral of p_k * d_phi_i / dx_j)
    B1_loc = -P' * (wq .* PX);
    B2_loc = -P' * (wq .* PY);

    % 7. Build D_loc (Integral of pressure basis functions)
    D_loc = P' * wq;

end
