clearvars; close all

meshRoot = fullfile(pwd, 'MESH');
resultDir = fullfile(pwd, 'results');

if ~exist(resultDir, 'dir')
    mkdir(resultDir);
end

% Mesh folders
meshTypes = {'tri', 'tri_p', 'quad', 'quad_p', 'poly','poly_hang','diam'};

% Methods definition
methods = struct( ...
    'name', {'vem','minitri'}, ...
    'fun',  {@stokes_vem_fun, @stokes_tri_fun}, ...
    'validMeshes', {{'diam'}} ...
);

for m = 1:numel(methods)
    method = methods(m);

    for mt = 1:numel(method.validMeshes)
        meshType = method.validMeshes{mt};
        meshPath = fullfile(meshRoot, meshType);

        mesh_files = dir(fullfile(meshPath, '*.txt'));
        if isempty(mesh_files)
            warning('No mesh files found in %s', meshPath);
            continue
        end

        results = cell(1, numel(mesh_files));

        fprintf('\n=== %s on %s meshes ===\n', method.name, meshType)

        for k = 1:numel(mesh_files)
            meshfile = fullfile(mesh_files(k).folder, mesh_files(k).name);
            fprintf('Processing %s...\n', meshfile)
            results{k} = method.fun(meshfile);
        end

        saveName = sprintf('%s_%s_result.mat', method.name, meshType);
        save(fullfile(resultDir, saveName), 'results');
    end
end
