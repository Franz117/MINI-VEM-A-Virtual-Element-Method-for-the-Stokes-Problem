clearvars; close all

resultDir = fullfile(pwd, 'results');
files = dir(fullfile(resultDir, '*_result.mat'));

for f = 1:numel(files)
    fprintf('\n\n=============================================\n');
    fprintf('Post-processing: %s\n', files(f).name);
    fprintf('=============================================\n');

    load(fullfile(files(f).folder, files(f).name), 'results');

    postprocess_results(results, files(f).name);
end
