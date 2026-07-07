clc; clear; close all;

% =====================================================================
% Plot CEC2017 convergence curves in the same style as the CEC2022 D10/D20
% figures.
% =====================================================================

%% ========================= CONFIGURATION =============================
result_mat = fullfile('results_cec2017', ...
    'CEC2017_BMJAYA_strong6_D50_N50_FEs500000_Runs30.mat');

if ~exist(result_mat, 'file')
    error(['Result MAT file not found: %s\n', ...
        'Run run_cec2017_bmjaya_strong6_curves.m first.'], result_mat);
end

show_iqr = false;
export_png = true;
export_pdf = true;
show_inset = true;
inset_fe_fraction = [0.75, 1.00];

%% ============================= LOAD ==================================
load(result_mat, 'config', 'curve_error_sampled');

alg_names = config.alg_names;
func_list = config.func_list;
fes_show = config.fes_show;

nFunc = numel(func_list);
nAlg = numel(alg_names);

titlestr = cell(1, nFunc);
panel_lab = cell(1, nFunc);
for p = 1:nFunc
    titlestr{p} = sprintf('F%d', func_list(p));
    panel_lab{p} = sprintf('(%c)', 'a' + p - 1);
end

marker_idx = 1:250:numel(fes_show);
[colors, line_styles] = fixed_algorithm_styles(alg_names);

bmjaya_idx = find(strcmpi(alg_names, 'BMJAYA'), 1);
if ~isempty(bmjaya_idx)
    plot_order = [setdiff(1:nAlg, bmjaya_idx, 'stable'), bmjaya_idx];
else
    plot_order = 1:nAlg;
end

%% ============================= PLOT ==================================
figure('Color','w','Position',[120 40 1100 1080]);
tl = tiledlayout(3, 2, 'Padding', 'loose', 'TileSpacing', 'loose');

for p = 1:nFunc
    ax = nexttile(tl, p);
    hold(ax, 'on');
    box(ax, 'on');
    grid(ax, 'on');

    med_curves = nan(nAlg, numel(fes_show));

    for order_id = 1:nAlg
        alg_id = plot_order(order_id);
        run_curves = double(squeeze(curve_error_sampled(:, p, alg_id, :)));
        med_curve = median(run_curves, 1);
        med_curves(alg_id, :) = med_curve;

        if show_iqr
            q1 = percentile_over_runs(run_curves, 0.25);
            q3 = percentile_over_runs(run_curves, 0.75);
            fill(ax, [fes_show, fliplr(fes_show)], ...
                [q1, fliplr(q3)], colors(alg_id, :), ...
                'FaceAlpha', 0.08, 'EdgeColor', 'none', ...
                'HandleVisibility', 'off');
        end

        semilogy(ax, fes_show, med_curve, line_styles{alg_id}, ...
            'LineWidth', 1.25, ...
            'MarkerIndices', marker_idx, ...
            'MarkerSize', 5, ...
            'Color', colors(alg_id, :), ...
            'DisplayName', alg_names{alg_id});
    end

    title(ax, titlestr{p});
    xlim(ax, [1, config.MaxFEs]);
    xlabel(ax, 'Function Evaluations (FEs)');
    ylabel(ax, 'Objective error |f(x)-f^*|');
    set(ax, 'FontSize', 9);

    text(ax, 0.50, -0.17, panel_lab{p}, ...
        'Units', 'normalized', ...
        'HorizontalAlignment', 'center', ...
        'FontWeight', 'bold');

    if show_inset && func_list(p) ~= 4
        add_late_stage_inset(ax, fes_show, med_curves, alg_names, colors, ...
            plot_order, config.MaxFEs, inset_fe_fraction);
    end

    legend(ax, 'Location', 'northeast', 'FontSize', 7);
end

if export_png
    exportgraphics(gcf, sprintf( ...
        'CEC2017_BMJAYA_strong6_convergence_D%d.png', config.D), ...
        'Resolution', 600);
end

if export_pdf
    exportgraphics(gcf, sprintf( ...
        'CEC2017_BMJAYA_strong6_convergence_D%d.pdf', config.D), ...
        'ContentType', 'vector');
end

%% ========================== LOCAL FUNCTIONS ==========================
function [colors, line_styles] = fixed_algorithm_styles(alg_names)
nAlg = numel(alg_names);
colors = zeros(nAlg, 3);
line_styles = cell(1, nAlg);

for alg_id = 1:nAlg
    name = upper(strtrim(alg_names{alg_id}));
    switch name
        case 'PSO'
            colors(alg_id, :) = [0.000, 0.447, 0.741];
            line_styles{alg_id} = '-o';
        case 'GWO'
            colors(alg_id, :) = [0.850, 0.325, 0.098];
            line_styles{alg_id} = '-^';
        case 'WSO'
            colors(alg_id, :) = [0.929, 0.694, 0.125];
            line_styles{alg_id} = '-x';
        case 'BMJAYA'
            colors(alg_id, :) = [0.494, 0.184, 0.556];
            line_styles{alg_id} = '-*';
        case 'DE'
            colors(alg_id, :) = [0.466, 0.674, 0.188];
            line_styles{alg_id} = '-d';
        case 'JAYA'
            colors(alg_id, :) = [0.301, 0.745, 0.933];
            line_styles{alg_id} = '-s';
        case 'SBOA'
            colors(alg_id, :) = [0.635, 0.078, 0.184];
            line_styles{alg_id} = '-v';
        case 'DOA'
            colors(alg_id, :) = [0.000, 0.200, 0.600];
            line_styles{alg_id} = '-p';
        otherwise
            fallback = lines(nAlg);
            colors(alg_id, :) = fallback(alg_id, :);
            line_styles{alg_id} = '-o';
    end
end
end

function add_late_stage_inset(parent_ax, fes_show, med_curves, alg_names, colors, ...
    plot_order, MaxFEs, inset_fe_fraction)
fig = ancestor(parent_ax, 'figure');
parent_pos = parent_ax.Position;
rel_pos = [0.30, 0.50, 0.36, 0.30];
inset_pos = [ ...
    parent_pos(1) + rel_pos(1) * parent_pos(3), ...
    parent_pos(2) + rel_pos(2) * parent_pos(4), ...
    rel_pos(3) * parent_pos(3), ...
    rel_pos(4) * parent_pos(4)];

inset_ax = axes(fig, 'Position', inset_pos);
hold(inset_ax, 'on');
box(inset_ax, 'on');
grid(inset_ax, 'on');

x_left = max(1, round(inset_fe_fraction(1) * MaxFEs));
x_right = max(x_left + 1, round(inset_fe_fraction(2) * MaxFEs));
idx = fes_show >= x_left & fes_show <= x_right;
if nnz(idx) < 2
    idx = fes_show >= max(1, round(0.5 * MaxFEs));
end

nAlg = size(med_curves, 1);
bmjaya_idx = find(strcmpi(alg_names, 'BMJAYA'), 1);
late_score = median(med_curves(:, idx), 2, 'omitnan');
[~, order] = sort(late_score, 'ascend');
focus_count = min(4, nAlg);
focus_idx = order(1:focus_count);
if ~isempty(bmjaya_idx) && ~ismember(bmjaya_idx, focus_idx)
    focus_idx(end) = bmjaya_idx;
end
inset_order = plot_order(ismember(plot_order, focus_idx));

for k = 1:numel(inset_order)
    alg_id = inset_order(k);
    if ~isempty(bmjaya_idx) && alg_id == bmjaya_idx
        inset_width = 2.00;
    else
        inset_width = 1.15;
    end

    plot(inset_ax, fes_show(idx), med_curves(alg_id, idx), ...
        '-', 'LineWidth', inset_width, 'Color', colors(alg_id, :));
end

y_local = med_curves(focus_idx, idx);
y_local = y_local(isfinite(y_local) & y_local > 0);
x_zoom_lim = [min(fes_show(idx)), max(fes_show(idx))];
y_zoom_lim = [];
if ~isempty(y_local)
    y_min = min(y_local);
    y_max = max(y_local);
    if y_max <= y_min
        y_max = y_min + max(1e-12, abs(y_min) * 0.1);
    end
    y_span = y_max - y_min;
    y_zoom_lim = [y_min - 0.10 * y_span, y_max + 0.10 * y_span];
    ylim(inset_ax, y_zoom_lim);
end
xlim(inset_ax, x_zoom_lim);
set(inset_ax, 'FontSize', 7, 'XTick', [], 'LineWidth', 0.7, 'Color', 'w');

if ~isempty(y_zoom_lim)
    y_box_lim = data_positive_box_limit(parent_ax, y_zoom_lim, y_local);
    draw_zoom_box(parent_ax, x_zoom_lim, y_box_lim);
    add_zoom_arrow(fig, parent_ax, inset_pos, x_zoom_lim, y_box_lim);
end
end

function y_box_lim = data_positive_box_limit(parent_ax, y_zoom_lim, y_local)
y_box_lim = y_zoom_lim;
if strcmpi(parent_ax.YScale, 'log')
    min_positive = min(y_local(y_local > 0));
    y_box_lim(1) = max(y_box_lim(1), min_positive * 0.80);
end
if y_box_lim(2) <= y_box_lim(1)
    y_box_lim(2) = y_box_lim(1) + max(1e-12, abs(y_box_lim(1)) * 0.1);
end
end

function draw_zoom_box(parent_ax, x_zoom_lim, y_box_lim)
x_box = [x_zoom_lim(1), x_zoom_lim(2), x_zoom_lim(2), x_zoom_lim(1), x_zoom_lim(1)];
y_box = [y_box_lim(1), y_box_lim(1), y_box_lim(2), y_box_lim(2), y_box_lim(1)];
plot(parent_ax, x_box, y_box, '--', ...
    'Color', [0.25, 0.25, 0.25], ...
    'LineWidth', 0.8, ...
    'HandleVisibility', 'off');
end

function add_zoom_arrow(fig, parent_ax, inset_pos, x_zoom_lim, y_box_lim)
arrow_start = [ ...
    inset_pos(1) + 0.96 * inset_pos(3), ...
    inset_pos(2) + 0.10 * inset_pos(4)];
target_x = x_zoom_lim(1) + 0.72 * (x_zoom_lim(2) - x_zoom_lim(1));
target_y = y_box_lim(1) + 0.50 * (y_box_lim(2) - y_box_lim(1));
arrow_end = data_to_figure(parent_ax, target_x, target_y);

annotation(fig, 'arrow', ...
    [arrow_start(1), arrow_end(1)], ...
    [arrow_start(2), arrow_end(2)], ...
    'Color', [0.25, 0.25, 0.25], ...
    'LineWidth', 0.8, ...
    'HeadLength', 6, ...
    'HeadWidth', 6);
end

function pt = data_to_figure(ax, x, y)
ax_pos = ax.Position;
x_lim = ax.XLim;
y_lim = ax.YLim;
if strcmpi(ax.XScale, 'log')
    x = log10(max(x, realmin));
    x_lim = log10(max(x_lim, realmin));
end
if strcmpi(ax.YScale, 'log')
    y = log10(max(y, realmin));
    y_lim = log10(max(y_lim, realmin));
end
x_frac = (x - x_lim(1)) / (x_lim(2) - x_lim(1));
y_frac = (y - y_lim(1)) / (y_lim(2) - y_lim(1));
x_frac = min(max(x_frac, 0), 1);
y_frac = min(max(y_frac, 0), 1);
pt = [ax_pos(1) + x_frac * ax_pos(3), ...
      ax_pos(2) + y_frac * ax_pos(4)];
end

function qv = percentile_over_runs(A, q)
A = sort(A, 1);
n = size(A, 1);
if n == 1
    qv = A;
    return;
end
pos = 1 + (n - 1) * q;
lo = floor(pos);
hi = ceil(pos);
w = pos - lo;
if lo == hi
    qv = A(lo, :);
else
    qv = (1 - w) .* A(lo, :) + w .* A(hi, :);
end
end
