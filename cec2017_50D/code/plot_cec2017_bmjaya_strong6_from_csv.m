clc; clear; close all;

% =====================================================================
% Plot CEC2017 final-performance comparison from CSV summary files.
%
% Data source:
%   E:\桌面\一修\结果\CEC2017_<ALG>_50D.csv
%
% The selected functions are the six CEC2017 functions where BMJAYA has
% the strongest ranks in the combined 50D MaxFEs=500000 result table.
% =====================================================================

%% ========================= CONFIGURATION =============================
data_dir = fileparts(mfilename('fullpath'));

alg_names = {'PSO','GWO','WSO','DE','JAYA','SBOA','DOA','BMJAYA'};
func_list = [3, 4, 12, 13, 18, 25];

D = 50;
MaxFEs = 500000;

export_png = true;
export_pdf = true;

[colors, markers] = fixed_algorithm_styles(alg_names);

%% ============================= LOAD ==================================
nAlg = numel(alg_names);
nFunc = numel(func_list);

mean_error = nan(nFunc, nAlg);
std_error = nan(nFunc, nAlg);
best_error = nan(nFunc, nAlg);

for alg_id = 1:nAlg
    alg = alg_names{alg_id};
    csv_path = fullfile(data_dir, sprintf('CEC2017_%s_50D.csv', alg));
    assert(isfile(csv_path), 'Cannot find CSV file: %s', csv_path);

    T = readtable(csv_path, 'TextType', 'string');

    for p = 1:nFunc
        func_id = func_list(p);
        row = T(strcmp(T.Function, sprintf('F%d', func_id)), :);
        assert(height(row) == 1, 'Cannot find F%d in %s.', func_id, csv_path);

        fopt = 100 * func_id;
        mean_error(p, alg_id) = abs(row.Mean - fopt) + 1e-12;
        std_error(p, alg_id) = row.Std;
        best_error(p, alg_id) = abs(row.Best - fopt) + 1e-12;
    end
end

%% ============================= PLOT ==================================
figure('Color', 'w', 'Position', [120 30 1100 1350]);
tl = tiledlayout(3, 2, 'Padding', 'loose', 'TileSpacing', 'loose');

for p = 1:nFunc
    ax = nexttile(tl, p);
    hold(ax, 'on');
    box(ax, 'on');
    grid(ax, 'on');

    x = 1:nAlg;

    for alg_id = 1:nAlg
        plot(ax, x(alg_id), mean_error(p, alg_id), ...
            markers{alg_id}, ...
            'LineStyle', 'none', ...
            'MarkerSize', 8, ...
            'LineWidth', 1.4, ...
            'Color', colors(alg_id, :), ...
            'MarkerFaceColor', marker_face_color(alg_names{alg_id}, colors(alg_id, :)), ...
            'DisplayName', alg_names{alg_id});
    end

    % Thin guide line for scanning across algorithms.
    plot(ax, x, mean_error(p, :), '-', ...
        'Color', [0.35, 0.35, 0.35], ...
        'LineWidth', 0.7, ...
        'HandleVisibility', 'off');

    set(ax, 'YScale', 'log');
    xlim(ax, [0.5, nAlg + 0.5]);
    set(ax, 'XTick', x, 'XTickLabel', alg_names, 'FontSize', 9);
    xtickangle(ax, 35);

    title(ax, sprintf('F%d', func_list(p)));
    ylabel(ax, 'Mean error |f(x)-f^*|');

    legend(ax, 'Location', 'northeast', 'FontSize', 7);

    text(ax, 0.04, 0.08, sprintf('(%c)', 'a' + p - 1), ...
        'Units', 'normalized', ...
        'HorizontalAlignment', 'left', ...
        'FontWeight', 'bold', ...
        'BackgroundColor', 'w', ...
        'Margin', 1);
end

if export_png
    exportgraphics(gcf, ...
        sprintf('CEC2017_BMJAYA_strong6_performance_D%d_MaxFEs%d.png', D, MaxFEs), ...
        'Resolution', 600);
end

if export_pdf
    exportgraphics(gcf, ...
        sprintf('CEC2017_BMJAYA_strong6_performance_D%d_MaxFEs%d.pdf', D, MaxFEs), ...
        'ContentType', 'vector');
end

%% ========================== LOCAL FUNCTIONS ==========================
function [colors, markers] = fixed_algorithm_styles(alg_names)
nAlg = numel(alg_names);
colors = zeros(nAlg, 3);
markers = cell(1, nAlg);

for alg_id = 1:nAlg
    name = upper(strtrim(alg_names{alg_id}));

    switch name
        case 'PSO'
            colors(alg_id, :) = [0.000, 0.447, 0.741];
            markers{alg_id} = 'o';
        case 'GWO'
            colors(alg_id, :) = [0.850, 0.325, 0.098];
            markers{alg_id} = '^';
        case 'WSO'
            colors(alg_id, :) = [0.929, 0.694, 0.125];
            markers{alg_id} = 'x';
        case 'DE'
            colors(alg_id, :) = [0.466, 0.674, 0.188];
            markers{alg_id} = 'd';
        case 'JAYA'
            colors(alg_id, :) = [0.301, 0.745, 0.933];
            markers{alg_id} = 's';
        case 'SBOA'
            colors(alg_id, :) = [0.635, 0.078, 0.184];
            markers{alg_id} = 'v';
        case 'DOA'
            colors(alg_id, :) = [0.000, 0.200, 0.600];
            markers{alg_id} = 'p';
        case 'BMJAYA'
            colors(alg_id, :) = [0.494, 0.184, 0.556];
            markers{alg_id} = '*';
        otherwise
            fallback = lines(nAlg);
            colors(alg_id, :) = fallback(alg_id, :);
            markers{alg_id} = 'o';
    end
end
end

function face = marker_face_color(alg_name, color)
if any(strcmpi(alg_name, {'PSO','GWO','DE','JAYA','SBOA','DOA','BMJAYA'}))
    face = color;
else
    face = 'none';
end
end
