function save_mesh(mesh, filename)
%SAVE_MESH  Writes mesh structure to file in required format

    fid = fopen(filename, 'w');
    
    if fid == -1
        error('Cannot open file.');
    end
    
    %% -------------------------------------------------
    % Domain type
    %% -------------------------------------------------
    fprintf(fid, '# domain type\n');
    fprintf(fid, '%s\n', mesh.type);
    
    %% -------------------------------------------------
    % Nodal coordinates
    %% -------------------------------------------------
    fprintf(fid, '# nodal coordinates: number of nodes followed by the coordinates\n');
    fprintf(fid, '%d\n', mesh.num_nodes);
    
    for i = 1:mesh.num_nodes
        fprintf(fid, '%.16f %.16f\n', ...
            mesh.coords(i,1), ...
            mesh.coords(i,2));
    end
    
    %% -------------------------------------------------
    % Element connectivity
    %% -------------------------------------------------
    fprintf(fid, '# element connectivity: number of elements followed by the connectivity\n');
    fprintf(fid, '%d\n', mesh.num_elements);
    
    for k = 1:mesh.num_elements
        nodes = mesh.connect{k};
        num_nodes_elem = length(nodes);
        
        % First print number of nodes in the element
        fprintf(fid, '%d ', num_nodes_elem);
        
        % Then print all node indices
        fprintf(fid, '%d ', nodes);
        
        fprintf(fid, '\n');
    end
    
    %% -------------------------------------------------
    % Boundary nodes
    %% -------------------------------------------------
    fprintf(fid, '# indices of nodes located on the boundary\n');
    
    for i = 1:length(mesh.boundary_nodes)
        fprintf(fid, '%d ', mesh.boundary_nodes(i));
    end
    fprintf(fid, '\n');
    
    fclose(fid);
end
