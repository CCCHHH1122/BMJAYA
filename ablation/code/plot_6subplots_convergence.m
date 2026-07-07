clc; clear; close all;
% ===== reproducibility =====


% ================== CONFIG ==================
func_list = [1, 2, 4, 6, 9, 11];
titlestr  = {'F1','F2','F4','F6','F9','F11'};
panel_lab = {'(a)','(b)','(c)','(d)','(e)','(f)'};

D = 20; 
N = 50;

MaxFEs  = 100000;   % real FE budget
MaxIter = 2000;     % shown iterations
lb = -100;
ub = 100;

K_req = 4;

% inset settings
inset_funcs = [1,2, 6, 9];     % 需要 inset 的函数（可改）
zoomL = 1800; 
zoomR = 2000;
show_mjaya_in_inset = true;   % inset 里是否画 MJAYA（true/false）

% 主图 marker 稀疏
iter_show = 0:MaxIter;
mk_main = 1:200:numel(iter_show);

% ================== FE -> Iter sampling (uniform, strictly increasing) ==================
fes_idx = round(linspace(1, MaxFEs, numel(iter_show)));
fes_idx(1) = 1;
fes_idx(end) = MaxFEs;
for i = 2:numel(fes_idx)
    if fes_idx(i) <= fes_idx(i-1)
        fes_idx(i) = min(MaxFEs, fes_idx(i-1) + 1);
    end
end

% ================== FIGURE LAYOUT ==================
fig = figure('Color','w','Position',[120 60 1000 940]);
tl = tiledlayout(3,2,'Padding','loose','TileSpacing','loose');



for p = 1:numel(func_list)

    func_id = func_list(p);   % ← 先定义 func_id

    baseSeed = 20260205;

    % ===== run algorithms (fixed & independent seeds) =====
    rng(baseSeed + func_id*100 + 1, 'twister');
    [~, c1] = alg_jaya_budget(func_id, D, N, MaxFEs, lb, ub);

    rng(baseSeed + func_id*100 + 2, 'twister');
    [~, c2] = alg_mjaya_fixedN_budget(func_id, D, N, MaxFEs, lb, ub, K_req);

    rng(baseSeed + func_id*100 + 3, 'twister');
    [~, c3] = alg_bmrjaya_budget(func_id, D, N, MaxFEs, lb, ub);

    rng(baseSeed + func_id*100 + 4, 'twister');
    [~, c4] = alg_bmjaya_fixedN_budget(func_id, D, N, MaxFEs, lb, ub, K_req);

    % ===== 后面绘图代码保持不变 =====
end


for p = 1:numel(func_list)
    func_id = func_list(p);

    % ===== run algorithms (must return curve of length MaxFEs) =====
    [~, c1] = alg_jaya_budget(func_id, D, N, MaxFEs, lb, ub);
    [~, c2] = alg_mjaya_fixedN_budget(func_id, D, N, MaxFEs, lb, ub, K_req);
    [~, c3] = alg_bmrjaya_budget(func_id, D, N, MaxFEs, lb, ub);
    [~, c4] = alg_bmjaya_fixedN_budget(func_id, D, N, MaxFEs, lb, ub, K_req);

    % map to Iter(0..2000)
    y1 = max(c1(fes_idx), 1e-12);
    y2 = max(c2(fes_idx), 1e-12);
    y3 = max(c3(fes_idx), 1e-12);
    y4 = max(c4(fes_idx), 1e-12);

    ax = nexttile(tl, p);
    hold(ax,'on'); box(ax,'on'); grid(ax,'on');

    % ===== axis style: log only when needed (rule-based, not by function ID) =====
    y_all = [y1(:);y2(:);y3(:);y4(:)];
    ymin_pos = min(y_all(y_all>0));
    dyn = max(y_all) / max(ymin_pos, 1e-12);
    use_log = dyn > 100;   % 超过2个数量级 -> log

    if use_log
        semilogy(ax, iter_show, y1, '-d','LineWidth',1.2,'MarkerIndices',mk_main);
        semilogy(ax, iter_show, y2, '-*','LineWidth',1.2,'MarkerIndices',mk_main);
        semilogy(ax, iter_show, y3, '-p','LineWidth',1.2,'MarkerIndices',mk_main);
        semilogy(ax, iter_show, y4, '-s','LineWidth',1.2,'MarkerIndices',mk_main);
    else
        plot(ax, iter_show, y1, '-d','LineWidth',1.2,'MarkerIndices',mk_main);
        plot(ax, iter_show, y2, '-*','LineWidth',1.2,'MarkerIndices',mk_main);
        plot(ax, iter_show, y3, '-p','LineWidth',1.2,'MarkerIndices',mk_main);
        plot(ax, iter_show, y4, '-s','LineWidth',1.2,'MarkerIndices',mk_main);

        % 线性图：聚焦平台段（避免初始大值压扁）
        idxP = iter_show >= 200;
        yP = [y1(idxP) y2(idxP) y3(idxP) y4(idxP)];
        yminP = min(yP(:)); ymaxP = max(yP(:));
        pad = 0.15*(ymaxP - yminP + eps);
        ylim(ax, [yminP - pad, ymaxP + pad]);
    end

    title(ax, titlestr{p});
    xlim(ax,[0 MaxIter]);
    set(ax,'XTick',0:400:2000);
    xlabel(ax,'Iteration');
    ylabel(ax,'Best score obtained so far');

    % legend only on right column
     if ismember(p,[1,2,3,4,5,6])
          legend(ax,'JAYA','MJAYA','BJAYA','BMJAYA','Location','northeast');
     end
 

    % ===== (a)-(f) at the bottom center (outside axes) =====
    text(ax, 0.50, -0.20, panel_lab{p}, 'Units','normalized', ...
        'HorizontalAlignment','center', 'VerticalAlignment','top', ...
        'FontSize',10, 'FontWeight','bold');

    % ================== inset (no markers) ==================
    if ismember(func_id, inset_funcs)
        idxZ = iter_show >= zoomL & iter_show <= zoomR;
        itz = iter_show(idxZ);
        z1 = y1(idxZ); z2 = y2(idxZ); z3 = y3(idxZ); z4 = y4(idxZ);

        % inset position relative to tile
        pos = ax.Position;
        % --- inset size & position (smaller, top-left to avoid legend) ---
        % --- inset size & position (right-middle, like the example figure) ---
        insetW = pos(3)*0.25;
        insetH = pos(4)*0.20;

        % horizontal: 靠右，但给 legend 留空间
        marginR = pos(3)*0.4;
        insetL  = pos(1) + pos(3) - insetW - marginR;

        % vertical: 中上部（不是顶到最上）
        centerY = pos(2) + pos(4)*0.55;
        insetB  = centerY - insetH/2;


        ax2 = axes('Position',[insetL insetB insetW insetH], 'Color','w');
        hold(ax2,'on'); box(ax2,'on'); grid(ax2,'on');
        % make inset non-interactive (avoid hit-test bug)
        set(ax2,'HitTest','off','PickableParts','none');
        set(ax2.Children,'HitTest','off','PickableParts','none');

        % inset: pure lines, no markers
        plot(ax2, itz, z1, '-', 'LineWidth',1.05);
        if show_mjaya_in_inset
            plot(ax2, itz, z2, '-', 'LineWidth',1.05);
        end
        plot(ax2, itz, z3, '-', 'LineWidth',1.05);
        plot(ax2, itz, z4, '-', 'LineWidth',1.25);  % BMJAYA emphasized

        set(ax2,'XLim',[zoomL zoomR]);
        set(ax2,'XTick',[]);
        set(ax2,'XTickLabel',[]);

        set(ax2,'FontSize',8);

        % inset y-range: focus on best three curves to avoid MJAYA blowing up scale
        z_best = [z1(:); z3(:); z4(:)];
        zmin = min(z_best); zmax = max(z_best);
        pad = 0.10*(zmax - zmin + eps);
        ylim(ax2, [zmin - pad, zmax + pad]);

        % fewer y ticks (clean)
        % fewer y ticks (clean) - robust & strictly increasing
%         yt = linspace(zmin, zmax, 3);
%         yt = unique(yt, 'stable');          % 去掉重复
%         if numel(yt) < 2                    % 极端：zmin == zmax
%           yt = [zmin, zmin + max(1e-6, 0.01*abs(zmin) + 1e-6)];
%         end
%         set(ax2,'YTick',yt);
        % ---- inset y-ticks: integers only (robust) ----
        yminI = floor(zmin);
        ymaxI = ceil(zmax);

        if yminI == ymaxI
          yt = [yminI, yminI+1];     % 防止只有一个刻度
        else
        % 最多 3 个整数刻度（和你原来一样简洁）
          yt = round(linspace(yminI, ymaxI, min(3, ymaxI-yminI+1)));
          yt = unique(yt, 'stable'); % 确保严格递增
        end

        set(ax2,'YTick', yt);
        set(ax2,'YTickLabel', string(yt));   % 明确显示为整数


        % connector line (thin, no arrow)
        xA = 1650;
        yA = interp1(iter_show, y4, xA, 'linear','extrap');
        xB = 2000;
        yB = interp1(iter_show, y4, 2000, 'linear','extrap');
        plot(ax, [xA xB], [yA yB], 'k-', 'LineWidth',0.8, 'HandleVisibility','off');
    end
end

% ================== DISABLE INTERACTIVITY (fix your version) ==================
% (1) disable on all axes (main + inset)
% axs = findall(gcf,'Type','axes');
% for k = 1:numel(axs)
%     try
%         disableDefaultInteractivity(axs(k));   % your MATLAB requires axes handle
%     catch
%     end
% end
% 
% % (2) turn off datatips + common mouse callbacks that trigger the bug
% try, datacursormode(gcf,'off'); catch, end
% try, set(gcf,'WindowButtonMotionFcn',[]); catch, end
% try, set(gcf,'WindowScrollWheelFcn',[]); catch, end
% try, set(gcf,'WindowButtonDownFcn',[]); catch, end

% ================== EXPORT (optional) ==================
% exportgraphics(gcf,'convergence_F1F3F5F7F9F11_D20.png','Resolution',300);
