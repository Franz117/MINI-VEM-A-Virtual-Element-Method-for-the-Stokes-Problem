function make_mesh_triang(n, filename,perturb_factor)

    % Triangulate 01 square with n triange per side
    [X, Y] = meshgrid(linspace(0,1,n), linspace(0,1,n));
    points = [X(:), Y(:)];

    if perturb_factor > 0
        points = perturb(points, perturb_factor);
    end

    DT = delaunayTriangulation(points);
    figure;
    triplot(DT);
    
    % Extract from mesh
    elements = DT.ConnectivityList;    
    boundaryEdges = freeBoundary(DT);
    boundaryNodes = unique(boundaryEdges(:));

    % Write mesh file
    fid = fopen(filename, 'w');
     
    fprintf(fid, "# domain type\n");
    fprintf(fid, "Triangular\n");
    fprintf(fid, "# nodal coordinates: number of nodes followed by the coordinates\n");
    fprintf(fid, "%d\n", size(points, 1));
    for i = 1:size(points, 1)
        fprintf(fid, "%.16f %.16f\n", points(i, 1), points(i, 2));
    end
    fprintf(fid, "# element connectivity: number of elements followed by the connectivity\n");
    fprintf(fid, "%d\n", size(elements, 1));
    for i = 1:size(elements, 1)
        fprintf(fid, "%d %d %d %d\n", 3, elements(i, 1), elements(i, 2), elements(i, 3));
    end
    fprintf(fid, "# indices of nodes located on the boundary\n");
    fprintf(fid, "%d ", boundaryNodes);
    fprintf(fid, "\n");
    
    fclose(fid);
    
    fprintf("Mesh file '%s' generated successfully.\n", filename);
    
end
