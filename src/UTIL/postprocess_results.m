function postprocess_results(results, label)

n = numel(results);

data(n) = struct('h', [], 'err0_p', [], 'err0_u', [], 'err1_u', []);

for k = 1:n
    res = results{k};
    data(k).h       = res.h;
    data(k).err0_p  = res.err0_p;
    data(k).err0_u  = res.err0_u;
    data(k).err1_u  = res.err1_u;
end

% ---- sort by descending h ----
[~, idx] = sort([data.h], 'descend');
data = data(idx);

% ---- extract arrays ----
h       = [data.h];
err0_p  = [data.err0_p];
err0_u  = [data.err0_u];
err1_u  = [data.err1_u];

% ---- compute EOCs ----
eoc_p  = log(err0_p(1:end-1) ./ err0_p(2:end)) ./ log(h(1:end-1) ./ h(2:end));
eoc_u0 = log(err0_u(1:end-1) ./ err0_u(2:end)) ./ log(h(1:end-1) ./ h(2:end));
eoc_u1 = log(err1_u(1:end-1) ./ err1_u(2:end)) ./ log(h(1:end-1) ./ h(2:end));

% ---- print table ----
fprintf('    h       ||p||_L2    EOC     ||u||_L2    EOC     ||u||_H1    EOC\n');
fprintf('--------------------------------------------------------------------------\n');

for i = 1:length(h)
    if i < length(h)
        fprintf('%7.4f   %8.2e  %5.2f   %8.2e  %5.2f   %8.2e  %5.2f\n', ...
            h(i), err0_p(i), eoc_p(i), err0_u(i), eoc_u0(i), err1_u(i), eoc_u1(i));
    else
        fprintf('%7.4f   %8.2e         %8.2e         %8.2e\n', ...
            h(i), err0_p(i), err0_u(i), err1_u(i));
    end
end

% ---- plot ----
figure('Name', label, 'NumberTitle', 'off');
loglog(h, err0_p, '-o', 'LineWidth', 2, 'DisplayName', '||p - p_h||_{L^2}');
hold on
loglog(h, err0_u, '-s', 'LineWidth', 2, 'DisplayName', '||u - u_h||_{L^2}');
loglog(h, err1_u, '-^', 'LineWidth', 2, 'DisplayName', '||u - u_h||_{H^1}');

% Reference slopes
h_ref = [min(h), max(h)];
loglog(h_ref, h_ref, '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 1.5, 'DisplayName', 'O(h)');
loglog(h_ref, h_ref.^2, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.5, 'DisplayName', 'O(h^2)');

grid on
xlabel('h (mesh size)');
ylabel('Error');
title(strrep(label, '_', '\_'));
legend('Location', 'southeast');

ylim([1e-4 2e1]);   % <-- fix y-axis range here

end
