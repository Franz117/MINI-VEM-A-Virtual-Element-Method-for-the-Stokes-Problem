function [err0_u_loc, err1_u_loc, err0_p_loc] = get_ERRloc_quad(coords,U1,U2,P,sol)

    u1 = sol.u1;
    u2 = sol.u2;
    p = sol.p;
    du1dx = sol.du1dx;
    du1dy = sol.du1dy;
    du2dx = sol.du2dx;
    du2dy = sol.du2dy;

    geom = get_geom(coords);
    area = geom.area;
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
    Ps = [p1(xq,yq),p2(xq,yq),p3(xq,yq),p4(xq,yq), b(xq,yq)];
    PX = [p1x(xq,yq),p2x(xq,yq),p3x(xq,yq),p4x(xq,yq),bx(xq,yq)];
    PY = [p1y(xq,yq),p2y(xq,yq),p3y(xq,yq),p4y(xq,yq),by(xq,yq)];

    % 3. err 0 (L2) of u
    u1_vals = u1(xq,yq);
    u2_vals = u2(xq,yq);
    u1h_vals = Ps * U1;
    u2h_vals = Ps * U2;

    err0_u1_loc = wq' * ((u1_vals-u1h_vals).^2);
    err0_u2_loc = wq' * ((u2_vals-u2h_vals).^2);
    err0_u_loc = err0_u1_loc + err0_u2_loc;

    % 4. err 0 (L2) of p
    p_vals = p(xq,yq);
    ph_vals = Ps(:,1:4) * P;

    err0_p_loc = wq' * ((p_vals-ph_vals).^2);

    % 5. err 1 (H1) of u 
    u1x_vals = du1dx(xq,yq);
    u1y_vals = du1dy(xq,yq);
    u2x_vals = du2dx(xq,yq);
    u2y_vals = du2dy(xq,yq);
    u1hx_vals = PX * U1;
    u1hy_vals = PY * U1;
    u2hx_vals = PX * U2;
    u2hy_vals = PY * U2;

    err1_u1_loc = wq' * ((u1x_vals - u1hx_vals).^2 + (u1y_vals - u1hy_vals).^2);
    err1_u2_loc = wq' * ((u2x_vals - u2hx_vals).^2 + (u2y_vals - u2hy_vals).^2);
    err1_u_loc = err1_u1_loc + err1_u2_loc;

return





