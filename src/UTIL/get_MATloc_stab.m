function [K_loc,B1_loc,B2_loc,D_loc,h] = get_MATloc_stab(coords)

    geom = get_geom(coords);
    bar = geom.bar;
    h = geom.h;
    area = geom.area;

    p1 = @(x,y) 1 + 0*x + 0*y;
    p2 = @(x,y) (x-bar(1))/h  + 0*y;
    p3 = @(x,y) 0*x + (y-bar(2))/h;

    nv = size(coords,1);

    Q = get_quadrature(1, coords);
    D0 =[p1(Q.Points(:,1),Q.Points(:,2)), p2(Q.Points(:,1),Q.Points(:,2)) , p3(Q.Points(:,1),Q.Points(:,2)) ];
    means = Q.Weights'*D0/area;

    % A_loc
    [P_NABLA,G_TILDE] = compute_P_NABLA(coords,geom);
    A_loc = P_NABLA'*G_TILDE*P_NABLA;
    A_loc = padarray(A_loc , [1 1], 0, 'post');

    % DOFI_loc
    means_phi = means*P_NABLA;
    I = eye(nv);
    D = [p1(coords(:,1),coords(:,2)) , p2(coords(:,1),coords(:,2)) , p3(coords(:,1),coords(:,2)) ];
    DOFI_loc = (I-D*P_NABLA)'*(I-D*P_NABLA);
    DOFI_loc = DOFI_loc + means_phi'*means_phi;
    DOFI_loc = padarray(DOFI_loc , [1 1], 0, 'post');
    DOFI_loc(nv+1,1:end-1)  = -means_phi;
    DOFI_loc(1:end-1,nv+1) = -means_phi';
    DOFI_loc(nv+1,nv+1) = 1;

    % making of K_loc
    K_loc = A_loc + DOFI_loc;

    % making of B1_loc
    B1_loc = zeros(nv+1,nv);
    B1_loc(nv+1,:) = area*P_NABLA(2,:)/h;
    B1_loc(nv+1,:) = B1_loc(nv+1,:);%+ sum(I-D*P_NABLA)*area/h;
    B1_loc(1:nv,:) = B1_loc(1:nv,:) + (I-D*P_NABLA)*area/h; %tentativo di stabilizz
    
    % making of B2_loc
    B2_loc = zeros(nv+1,nv);
    B2_loc(nv+1,:) = area*P_NABLA(3,:)/h;
    B2_loc(nv+1,:) = B2_loc(nv+1,:);%+ sum(I-D*P_NABLA)*area/h;
    B2_loc(1:nv,:) = B2_loc(1:nv,:) + (I-D*P_NABLA)*area/h; %tentativo di stabilizz

    % making of B1_bord_loc and B2_bord_loc
    Q = get_border_quadrature_AI(geom, coords);
    D1 =[p1(Q.Points(:,1),Q.Points(:,2)), p2(Q.Points(:,1),Q.Points(:,2)) , p3(Q.Points(:,1),Q.Points(:,2)) ];
    Psi = D1 * P_NABLA;
    nq = numel(Q.Points)/2;
    Phi = get_PHI_GAUSS_2(nq,nv);
    %Phi = Psi;%stabilizzazione (questa forma si può calcolare esattamente)
    W1 = spdiags(Q.Weights1,0,length(Q.Weights1),length(Q.Weights1)); 
    W2 = spdiags(Q.Weights2,0,length(Q.Weights2),length(Q.Weights2)); 
    B1_bord = Phi' * W1 * Psi; 
    B2_bord = Phi' * W2 * Psi; 
    B1_loc(1:nv,1:nv) = B1_loc(1:nv,1:nv) - B1_bord;
    B2_loc(1:nv,1:nv) = B2_loc(1:nv,1:nv) - B2_bord;

    % making of D_loc
    D_loc = area*means_phi;

return

function Phi = get_PHI_GAUSS_2(nq,nv)
% calculate by hand matrix where each colum is the basis function evaulated
% in the noed of a gauss 2 quadrature; very ugly.

    Phi = zeros(nq,nv);
    gv1 = (3-sqrt(3))/6;
    gv2 = (3+sqrt(3))/6;
    Phi(nq-1,1)=gv1;
    Phi(nq,1)=gv2;
    Phi(1,1)=gv2;
    Phi(2,1)=gv1;
    for i=2:nv
        Phi(i*2-3,i)=gv1;
        Phi(i*2-2,i)=gv2;
        Phi(i*2-1,i)=gv2;
        Phi(i*2,i)=gv1;
    end

return