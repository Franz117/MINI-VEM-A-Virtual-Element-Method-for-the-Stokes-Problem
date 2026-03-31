function [K_loc,B1_loc,B2_loc,D_loc,h] = get_MATloc_quad(coords)

    geom = get_geom(coords);
    h = geom.h;
    area = geom.area;
    edges = geom.edges;
    xv = coords(:,1);
    yv = coords(:,2);

    x1 = coords(1,1); y1 = coords(1,2);
    x2 = coords(2,1); y2 = coords(2,2);
    x3 = coords(3,1); y3 = coords(3,2);
    x4 = coords(4,1); y4 = coords(4,2);

    % 1. Basis functions
    J  = (x2-x1)*(y4-y1) - (y2-y1)*(x4-x1);
    s = @(x,y) ( (x-x1)*(y4-y1) - (y-y1)*(x4-x1) ) / J;
    t = @(x,y) (-(x-x1)*(y2-y1) + (y-y1)*(x2-x1) ) / J;
    sx =  (y4 - y1) / J;
    sy = -(x4 - x1) / J;
    tx = -(y2 - y1) / J;
    ty =  (x2 - x1) / J;

    p1 = @(x,y) (1 - s(x,y)).*(1 - t(x,y));
    p2 = @(x,y)      s(x,y) .*(1 - t(x,y));
    p3 = @(x,y)      s(x,y) .*      t(x,y);
    p4 = @(x,y) (1 - s(x,y)).*      t(x,y);
    %b = @(x,y) 16 *s(x,y).*(1 - s(x,y)).*t(x,y).*(1 - t(x,y));
    b = @(x,y) 8 * (1 + s(x,y) + t(x,y)) .* s(x,y).*(1-s(x,y)) .* t(x,y).*(1-t(x,y)); %BISHNU P. LAMICHHANE

    p1x = @(x,y) -(1-t(x,y))*sx - (1-s(x,y))*tx;
    p1y = @(x,y) -(1-t(x,y))*sy - (1-s(x,y))*ty;
    p2x = @(x,y)  (1-t(x,y))*sx - s(x,y)*tx;
    p2y = @(x,y)  (1-t(x,y))*sy - s(x,y)*ty;
    p3x = @(x,y)  t(x,y)*sx + s(x,y)*tx;
    p3y = @(x,y)  t(x,y)*sy + s(x,y)*ty;
    p4x = @(x,y) -t(x,y)*sx + (1-s(x,y))*tx;
    p4y = @(x,y) -t(x,y)*sy + (1-s(x,y))*ty;
%     b1x = @(x,y) 16*((1-2*s(x,y)).*t(x,y).*(1-t(x,y))*sx ...
%                    + (1-2*t(x,y)).*s(x,y).*(1-s(x,y))*tx);
%     b1y = @(x,y) 16*((1-2*s(x,y)).*t(x,y).*(1-t(x,y))*sy ...
%                    + (1-2*t(x,y)).*s(x,y).*(1-s(x,y))*ty);
    b_s = @(s,t) 8 * ( t.*(1-t) ) .* ( s.*(1-s) + (1 + s + t).*(1 - 2*s) );
    b_t = @(s,t) 8 * ( s.*(1-s) ) .* ( t.*(1-t) + (1 + s + t).*(1 - 2*t) );
    bx = @(x,y) b_s(s(x,y), t(x,y)) * sx + b_t(s(x,y), t(x,y)) * tx;
    by = @(x,y) b_s(s(x,y), t(x,y)) * sy + b_t(s(x,y), t(x,y)) * ty;
    
    % 3. Quadrature (Weights include area)
    Q  = get_quadrature_quad(4, coords);
    xq = Q.Points(:,1);
    yq = Q.Points(:,2);
    wq = Q.Weights;
    P = [p1(xq,yq),p2(xq,yq),p3(xq,yq),p4(xq,yq)];
    PX = [p1x(xq,yq),p2x(xq,yq),p3x(xq,yq),p4x(xq,yq),bx(xq,yq)];
    PY = [p1y(xq,yq),p2y(xq,yq),p3y(xq,yq),p4y(xq,yq),by(xq,yq)];

    % 4. Build K_loc (Stiffness matrix: integral of grad_phi_i . grad_phi_j)
    K_loc = PX' * (wq .* PX) + PY' * (wq .* PY);

    % 5. Build B1_loc, B2_loc (Divergence: integral of p_k * d_phi_i / dx_j)
    B1_loc = -P' * (wq .* PX);
    B2_loc = -P' * (wq .* PY);

    % 7. Build D_loc (Integral of pressure basis functions)
    D_loc = P' * wq;

end
