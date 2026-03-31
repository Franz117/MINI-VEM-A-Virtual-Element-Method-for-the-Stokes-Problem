clearvars, close all

nn =[4,6,9,13,18,25,36,51,72];

for i=1:numel(nn)
    n=nn(i);
    make_mesh_quad1(n,1/(n*2.5));
end
