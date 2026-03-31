function means = compute_mean(coords, geom)
    % coords: Nx2 array of vertices [x y] in order

    bar = geom.bar;
    area = geom.area;
    h = geom.h;

    % Step 3: create constrained edges
    edges = [(1:size(coords,1))', [2:size(coords,1) 1]'];

    % Step 4: constrained Delaunay triangulation
    DT = delaunayTriangulation(coords, edges);
    inside = isInterior(DT);  % Keep only triangles inside the polygon
    TR = triangulation(DT.ConnectivityList(inside,:), DT.Points);

    % Step 5: initialize integrals
    int_p1 = 0;
    int_p2 = 0;
    int_p3 = 0;

    % Step 6: loop over triangles
    for k = 1:size(TR.ConnectivityList,1)
        % Get triangle vertices
        tri_indices = TR.ConnectivityList(k,:);
        tri_coords = TR.Points(tri_indices,:);

        % Area of triangle
        tri_area = polyarea(tri_coords(:,1), tri_coords(:,2));

        % Centroid of triangle
        xm = mean(tri_coords(:,1));
        ym = mean(tri_coords(:,2));

        % Evaluate functions at centroid
        p1_val = 1;
        p2_val = (xm - bar(1)) / h;
        p3_val = (ym - bar(2)) / h;

        % Add contribution
        int_p1 = int_p1 + p1_val * tri_area;
        int_p2 = int_p2 + p2_val * tri_area;
        int_p3 = int_p3 + p3_val * tri_area;
    end

    % Normalize by area
    int_p1 = int_p1 / area;
    int_p2 = int_p2 / area;
    int_p3 = int_p3 / area;
    means = [int_p1, int_p2, int_p3];
end
