function [err0_u_loc, err1_u_loc, err0_p_loc] = get_ERRloc_vem(coords,dof,U1,U2,P,sol)

    u1 = sol.u1;
    u2 = sol.u2;
    p = sol.p;
    du1dx = sol.du1dx;
    du1dy = sol.du1dy;
    du2dx = sol.du2dx;
    du2dy = sol.du2dy;

    geom = get_geom(coords);
    bar = geom.bar;
    h = geom.h;

    p1 = @(x,y) 1 + 0*x + 0*y;
    p2 = @(x,y) (x-bar(1))/h  + 0*y;
    p3 = @(x,y) 0*x + (y-bar(2))/h;

    Q = get_quadrature(8, coords);
    Xq = Q.Points(:,1);
    Yq = Q.Points(:,2);

    [P_NABLA,~] = compute_P_NABLA(coords,geom);
    D0 =[p1(Xq,Yq), p2(Xq,Yq) , p3(Xq,Yq) ];

    Psi = D0*P_NABLA;

    u1_vals = u1(Xq,Yq);
    u1h_vals = Psi * U1(dof);
    err0_u1_loc = sum(Q.Weights .* ((u1_vals - u1h_vals).^2));

    u2_vals = u2(Xq,Yq);
    u2h_vals = Psi * U2(dof);
    err0_u2_loc = sum(Q.Weights .* ((u2_vals - u2h_vals).^2));

    err0_u_loc = err0_u1_loc + err0_u2_loc;

    p_vals = p(Xq,Yq);
    ph_vals = Psi * P(dof);
    err0_p_loc = sum(Q.Weights .* ((p_vals - ph_vals).^2));

    du1dx_vals = du1dx(Xq,Yq);
    du1dy_vals = du1dy(Xq,Yq);
    grad_u1h_x = (P_NABLA(2,:) * U1(dof))/h;  % scalar
    grad_u1h_y = (P_NABLA(3,:) * U1(dof))/h;  % scalar
    err_x = du1dx_vals - grad_u1h_x;
    err_y = du1dy_vals - grad_u1h_y;
    err1_u1_loc = sum(Q.Weights .* (err_x.^2 + err_y.^2));

    du2dx_vals = du2dx(Xq,Yq);
    du2dy_vals = du2dy(Xq,Yq);
    grad_u2h_x = (P_NABLA(2,:) * U2(dof))/h;  % scalar
    grad_u2h_y = (P_NABLA(3,:) * U2(dof))/h;  % scalar
    err_x = du2dx_vals - grad_u2h_x;
    err_y = du2dy_vals - grad_u2h_y;
    err1_u2_loc = sum(Q.Weights .* (err_x.^2 + err_y.^2));
    
    err1_u_loc = err1_u1_loc + err1_u2_loc;


return





