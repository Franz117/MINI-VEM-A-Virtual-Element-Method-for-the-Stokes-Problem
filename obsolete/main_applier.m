% MIXED STOKES MULTIMETHOD APPLIER

clearvars; close all

%To Bulk apply

meshPath = fullfile(pwd, 'MESH', 'poly');
mesh_files  = dir(fullfile(meshPath, '*.txt'));

results = {};

for k = 1:length(mesh_files)
    meshfile = fullfile(mesh_files(k).folder, mesh_files(k).name);
    fprintf("Processing %s...\n", meshfile)
    results{k} = stokes_unvem_fun(meshfile);
end

%save("results/mini_tri_result.mat", "results");
%save("results//mini_tri_p_result.mat", "results");
%save("results/vem_tri_result.mat", "results");
%save("results/vem_tri_p_result.mat", "results");
%save("results/unvem_tri_result.mat", "results");
%save("results/unvem_tri_p_result.mat", "results");

%save("results/mini_quad_result.mat", "results");
%save("results//mini_quad_p_result.mat", "results");
%save("results/vem_quad_result.mat", "results");
%save("results/vem_quad_p_result.mat", "results");
%save("results/unvem_quad_result.mat", "results");
%save("results/unvem_quad_p_result.mat", "results");

%save("results/vem_poly_result.mat", "results");
save("results/unvem_poly_result.mat", "results");


% % Load the results
% load("results/mini_tri_result.mat", "results");
% n = numel(results);
% data(n) = struct('h', [], 'err0_p', [], 'err0_u', [], 'err1_u', []);
% for k = 1:n
%     res = results{k};
%     data(k).h       = res.h;
%     data(k).err0_p  = res.err0_p;
%     data(k).err0_u  = res.err0_u;
%     data(k).err1_u  = res.err1_u;
% end
% [~, idx] = sort([data.h]);
% data_sorted = data(idx);
% 
% % Extract data from data_sorted
% h       = [data_sorted.h];
% err0_p  = [data_sorted.err0_p];
% err0_u  = [data_sorted.err0_u];
% err1_u  = [data_sorted.err1_u];
% 
% 
% % Compute EOCs
% eoc_p  = log(err0_p(1:end-1) ./ err0_p(2:end)) ./ log(h(1:end-1) ./ h(2:end));
% eoc_u0 = log(err0_u(1:end-1) ./ err0_u(2:end)) ./ log(h(1:end-1) ./ h(2:end));
% eoc_u1 = log(err1_u(1:end-1) ./ err1_u(2:end)) ./ log(h(1:end-1) ./ h(2:end));
% 
% % Display in table format
% fprintf('    h       ||p||_L2    EOC     ||u||_L2    EOC     ||u||_H1    EOC\n');
% fprintf('--------------------------------------------------------------------------\n');
% for i = 1:length(h)
%     if i < length(h)
%         fprintf('%7.4f   %8.2e  %5.2f   %8.2e  %5.2f   %8.2e  %5.2f\n', ...
%             h(i), err0_p(i), eoc_p(i), err0_u(i), eoc_u0(i), err1_u(i), eoc_u1(i));
%     else
%         fprintf('%7.4f   %8.2e         %8.2e         %8.2e\n', ...
%             h(i), err0_p(i), err0_u(i), err1_u(i));
%     end
% end
% 
% % Plot errors
% figure;
% loglog(h, err0_p, '-o', 'LineWidth', 2, 'DisplayName', '||p - p_h||_{L^2}');
% hold on;
% loglog(h, err0_u, '-s', 'LineWidth', 2, 'DisplayName', '||u - u_h||_{L^2}');
% loglog(h, err1_u, '-^', 'LineWidth', 2, 'DisplayName', '||u - u_h||_{H^1}');
% 
% % Reference lines: O(h), O(h^2)
% h_ref = [min(h), max(h)];
% ref1 = h_ref;       % O(h)
% ref2 = h_ref.^2;    % O(h^2)
% 
% loglog(h_ref, ref1, '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 1.5, 'DisplayName', 'O(h)');
% loglog(h_ref, ref2, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.5, 'DisplayName', 'O(h^2)');
% 
% % Labels and formatting
% grid on;
% xlabel('h (mesh size)');
% ylabel('Error');
% title('Convergence of Mixed VEM Method');
% legend('Location', 'southeast');



