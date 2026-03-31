function [F1_loc,F2_loc] = get_Floc_vem(coords,f1,f2)

    nv = size(coords,1);
    geom = get_geom(coords);
    bar = geom.bar;
    area = geom.area;

    F1_loc = zeros(1,nv);
    F1_loc = [F1_loc, area*f1(bar(1),bar(2))];
    F2_loc = zeros(1,nv);
    F2_loc = [F2_loc, area*f2(bar(1),bar(2))];

return