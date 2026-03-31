function [K_loc,B1_loc,B2_loc,D_loc,h,F1_loc,F2_loc] = get_MATloc_unvem(coords,f1,f2)

    geom = get_geom(coords);
    h = geom.h;
   
    % 3. Quadrature (Weights include area)
    Q  = get_quadrature(6, coords);
    xq = Q.Points(:,1);
    yq = Q.Points(:,2);
    wq = Q.Weights;
    F1 = [f1(xq,yq)];
    F2 = [f2(xq,yq)];
    [P,PX,PY] = get_base_unvem(coords,xq,yq);
    P1 = P(:,1:end-1);

    % 4. Build K_loc (Stiffness matrix: integral of grad_phi_i . grad_phi_j)
    K_loc = PX' * (wq .* PX) + PY' * (wq .* PY);

    % 5. Build B1_loc, B2_loc (Divergence: integral of p_k * d_phi_i / dx_j)
    B1_loc = -P1' * (wq .* PX);
    B2_loc = -P1' * (wq .* PY);

    % 7. Build D_loc (Integral of pressure basis functions)
    D_loc = P1' * wq;

    % 8. F1_loc F2_loc
    F1_loc = P' * (wq .* F1);
    F2_loc = P' * (wq .* F2);

end
