function mesh = remove_nodes(mesh, remove_ids)
%REMOVE_NODES_FROM_MESH Removes specified nodes and updates mesh
%
% mesh        : input mesh structure
% remove_ids  : vector of node indices to remove
%
% returns updated mesh

    % --- sanity ---
    remove_ids = unique(remove_ids(:)); % ensure column & unique
    N = mesh.num_nodes;

    if any(remove_ids < 1) || any(remove_ids > N)
        error('remove_ids contains invalid node indices');
    end

    % --- mark nodes to remove ---
    to_remove = false(N,1);
    to_remove(remove_ids) = true;

    % --- keep nodes ---
    keep_ids = find(~to_remove);

    % --- build reindex map ---
    new_index = zeros(N,1);
    new_index(keep_ids) = 1:length(keep_ids);

    % --- update coordinates ---
    mesh.coords = mesh.coords(keep_ids,:);
    mesh.num_nodes = size(mesh.coords,1);

    % --- update connectivity ---
    new_connect = cell(size(mesh.connect));

    for e = 1:mesh.num_elements
        conn = mesh.connect{e};

        % remove deleted nodes
        conn = conn(~to_remove(conn));

        % reindex
        conn = new_index(conn);

        % optional: remove degenerate elements
        if numel(conn) < 3
            new_connect{e} = []; % mark for deletion
        else
            % optional: remove duplicate consecutive nodes
            conn = conn([true; diff(conn(:))~=0]);
            new_connect{e} = conn;
        end
    end

    % --- remove empty elements ---
    valid = ~cellfun(@isempty, new_connect);
    mesh.connect = new_connect(valid);
    mesh.num_elements = numel(mesh.connect);

    % --- update boundary nodes ---
    bnd = mesh.boundary_nodes;
    bnd = bnd(~to_remove(bnd));
    mesh.boundary_nodes = new_index(bnd);

end