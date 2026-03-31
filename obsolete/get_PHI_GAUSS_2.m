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