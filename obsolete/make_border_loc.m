function [B1_loc,B2_loc] = make_border_loc(coords)

    geom = get_geom(coords);
    bar = geom.bar;
    h = geom.h;

    p1 = @(x,y) 1 + 0*x + 0*y;
    p2 = @(x,y) (x-bar(1))/h  + 0*y;
    p3 = @(x,y) 0*x + (y-bar(2))/h;

    nv = size(coords,1);

    [P_NABLA,~] = compute_P_NABLA(coords,geom);

    B1_loc = zeros(nv+1,nv);
    B2_loc = zeros(nv+1,nv);

    %border correction with Cavalieri simpson 
    Q = get_border_quadrature(geom,coords);

    D1 =[p1(Q.Points(:,1),Q.Points(:,2)), p2(Q.Points(:,1),Q.Points(:,2)) , p3(Q.Points(:,1),Q.Points(:,2)) ];

    Phi_mat = zeros(nv, 2*nv);
    idx = 1:nv;
    Phi_mat(sub2ind(size(Phi_mat), idx, 2*idx - 1)) = 1;
    Phi_mat(sub2ind(size(Phi_mat), idx, 2*idx))= 0.5;
    non_first = idx > 1;
    Phi_mat(sub2ind(size(Phi_mat), idx(non_first), 2*idx(non_first)-2)) = 0.5;
    Phi_mat(1, end) = 0.5;
    
    Psi = D1 * P_NABLA;
    B1_bord = Phi_mat * (Q.Weights1 .* Psi);
    B2_bord = Phi_mat * (Q.Weights2 .* Psi);

    B1_loc(1:nv,1:nv) = B1_loc(1:nv,1:nv) -B1_bord;
    B2_loc(1:nv,1:nv) = B2_loc(1:nv,1:nv) -B2_bord;


return