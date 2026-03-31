function [err0_u_loc, err1_u_loc, err0_p_loc] = get_error_loc_proj(coords,dof,U1,U2,P,sol)

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
    area = geom.area;

    p1 = @(x,y) 1 + 0*x + 0*y;
    p2 = @(x,y) (x-bar(1))/h  + 0*y;
    p3 = @(x,y) 0*x + (y-bar(2))/h;

    Q = get_quadrature(2, coords);
    Xq = Q.Points(:,1);
    Yq = Q.Points(:,2);
    [P_NABLA,~] = compute_P_NABLA(coords,geom);

    p_diff = P(dof)-p(coords(:,1),coords(:,2));
    u1_diff = U1(dof)-u1(coords(:,1),coords(:,2));
    u2_diff = U2(dof)-u2(coords(:,1),coords(:,2));


    D0 =[p1(Xq,Yq), p2(Xq,Yq) , p3(Xq,Yq) ];
    Psi = D0*P_NABLA;

    err0_p_loc = Q.Weights'*sum(p_diff'.*Psi,2).^2;
    err0_u1_loc = Q.Weights'*sum(u1_diff'.*Psi,2).^2;
    err0_u2_loc = Q.Weights'*sum(u2_diff'.*Psi,2).^2;
    err0_u_loc = err0_u1_loc + err0_u2_loc;

    err1_u1_loc = norm(sum(u1_diff'.*P_NABLA(2:3,:),2)/h)*area;
    err1_u2_loc = norm(sum(u2_diff'.*P_NABLA(2:3,:),2)/h)*area;
    err1_u_loc = err1_u1_loc + err1_u2_loc;


return





