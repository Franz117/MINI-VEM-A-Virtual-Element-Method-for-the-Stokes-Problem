function mesh_maker(mesh_type, n, perturb_factor)
    if nargin < 3 || isequal(perturb_factor, "auto")
        perturb_factor = NaN; % placeholder
    end

    % Compute grid size
    switch mesh_type
        case "tri"
            n_grid = round(sqrt(n / 2)) + 1;
            actual_n = 2 * (n_grid - 1)^2;
            if isnan(perturb_factor)
                perturb_factor = 1 / (3 * (n_grid - 1));
            end

        case "quad"
            n_grid = round(sqrt(n)) + 1;
            actual_n = (n_grid - 1)^2;
            if isnan(perturb_factor)
                perturb_factor = 1 / (3 * (n_grid - 1));
            end

        case "poly"
            actual_n = n;
            n_grid = n; % just for filename
            perturb_factor = 0; % unused
        otherwise
            error("Unknown mesh type: %s", mesh_type);
    end

    suffix = "";
    if perturb_factor > 0 && (mesh_type == "tri" || mesh_type == "quad")
        suffix = "_p";
    end

    filename = sprintf("mesh_%s%s_%d.txt", mesh_type, suffix, actual_n);

    switch mesh_type
        case "tri"
            make_mesh_triang(n_grid, filename, perturb_factor);
        case "quad"
            make_mesh_quad(n_grid - 1, filename, perturb_factor);
        case "poly"
            make_mesh_poly(n, filename);
    end
end
