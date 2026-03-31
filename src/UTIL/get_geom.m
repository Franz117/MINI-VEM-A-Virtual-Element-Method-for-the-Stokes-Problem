function geom = get_geom(coords)

    %barycenter
    bar = mean(coords, 1);

    %element diamter
    h = max(pdist(coords, 'euclidean'));

    %edges lenghts
    coords = [coords; coords(1,:)]; %looping coords
    edges = sqrt(sum(diff(coords).^2, 2));
    border = sum(edges);

    %area
    area = 0.5 * abs(sum(coords(1:end-1,1) .* coords(2:end,2) - coords(2:end,1) .* coords(1:end-1,2)));

    %normals
    edges_vector = diff(coords); % Edge vectors (dx, dy)
    normals = [-edges_vector(:,2), edges_vector(:,1)]; % Rotate 90°
    normals = -normals ./ edges; % Normalize to unit vectors

    geom.bar = bar;
    geom.h = h;
    geom.edges = edges;
    geom.area = area;
    geom.border = border;
    geom.normals = normals;

end