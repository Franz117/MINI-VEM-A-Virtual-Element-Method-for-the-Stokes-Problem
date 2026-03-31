function [err0_u_loc, err1_u_loc, err0_p_loc] = get_ERRloc_tri(coords,U1,U2,P,sol)

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
    Ps = [p1(xq,yq),p2(xq,yq),p3(xq,yq), b(xq,yq)];
    PX = [p1x(xq,yq),p2x(xq,yq),p3x(xq,yq),bx(xq,yq)];
    PY = [p1y(xq,yq),p2y(xq,yq),p3y(xq,yq),by(xq,yq)];

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
    ph_vals = Ps(:,1:3) * P;

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





