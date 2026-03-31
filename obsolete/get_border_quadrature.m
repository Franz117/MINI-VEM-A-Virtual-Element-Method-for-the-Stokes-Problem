function Q = get_border_quadrature(geom,coords)

    edges = geom.edges;
    normals = geom.normals;

    midpoints = (coords(1:end-1, :) + coords(2:end, :)) / 2;
    last_midpoint = (coords(end, :) + coords(1, :)) / 2;
    all_midpoints = [midpoints; last_midpoint];
    new_coords = zeros(size(coords,1)*2, size(coords,2));
    new_coords(1:2:end, :) = coords;
    new_coords(2:2:end, :) = all_midpoints;
    Q.Points = new_coords;

    w1 = edges .* normals(:,1)/6; 
    w1(2:end) = w1(2:end) + w1(1:end-1);
    w1(1) = w1(1) + edges(end) * normals(end,1);
    insert_vals = edges .* normals(:,1) * (2/3);
    w1_interleaved = zeros(2 * length(w1), 1);
    w1_interleaved(1:2:end) = w1;
    w1_interleaved(2:2:end) = insert_vals;
    w1 = w1_interleaved;
    Q.Weights1 = w1;

    w2 = edges .* normals(:,2)/6; 
    w2(2:end) = w2(2:end) + w2(1:end-1);
    w2(1) = w2(1) + edges(end) * normals(end,2);
    insert_vals = edges .* normals(:,2) * (2/3);
    w2_interleaved = zeros(2 * length(w2), 1);
    w2_interleaved(1:2:end) = w2;
    w2_interleaved(2:2:end) = insert_vals;
    w2 = w2_interleaved;
    Q.Weights2 = w2;



return
