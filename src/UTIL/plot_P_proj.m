function plot_P_proj(mesh, P, title_str, use_surface)

ne = mesh.num_elements;
nv = mesh.num_nodes;

P_proj = zeros(size(P));
areas = zeros(nv,1);

for ie=1:ne 
    dof = mesh.connect{ie};
    coords = mesh.coords(dof,:);
    geom = get_geom(coords);
    area = geom.area;
    bar = geom.bar;
    h = geom.h;
    Ps = P(dof);

    p1 = @(x,y) 1 + 0*x + 0*y;
    p2 = @(x,y) (x-bar(1))/h  + 0*y;
    p3 = @(x,y) 0*x + (y-bar(2))/h;

    D = [p1(coords(:,1),coords(:,2)) , p2(coords(:,1),coords(:,2)) , p3(coords(:,1),coords(:,2)) ];
    [P_NABLA,~] = compute_P_NABLA(coords,geom);

    P_proj(dof) = P_proj(dof) + (Ps'*D*P_NABLA)'*area;
    areas(dof) = areas(dof) + area;
end

P_proj = P_proj./areas;
plot_P(mesh, P_proj, title_str, use_surface);


end
