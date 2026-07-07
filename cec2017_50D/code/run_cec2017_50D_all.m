function run_cec2017_50D_all()
clc; clear; close all;

% =====================================================================
% CEC2017 50D strict-MaxFEs benchmark runner for BMJAYA paper experiments.
%
% Active algorithms:
%   PSO, GWO, WSO, BMJAYA, DE, JAYA, DOA, SBOA
%
% Outputs:
%   - one CSV per algorithm
%   - combined detailed Rank CSV
%   - scientific-notation Rank CSV
%   - final SumRank / FinalRank CSV
%   - MAT file with sampled convergence curves
%
% Required official benchmark files:
%   cec17_func.mexw64
%   input_data\
% =====================================================================

%% ========================= USER CONFIG ===============================
CEC_DIR = 'E:\桌面\project\cec2017';
OUT_DIR = fullfile(pwd, 'results_cec2017');

Dims = [50];
N = 50;
MaxFEs = 5000000;
Runs = 30;

% Standard paper-ready CEC2017 set:
% F2 is commonly excluded from the CEC2017 bound-constrained comparison
% because of known instability in the original benchmark implementation.
FuncIDs = [1, 3:30];

% Representative functions for convergence figures.
CurveFuncIDs = [1, 3, 10, 15, 20, 30];

Algorithms = { ...
    %struct('name','PSO'   , 'fn', @alg_pso_weak_cec17); ...
    %struct('name','GWO'   , 'fn', @alg_gwo_cec17); ...
    struct('name','WSO'   , 'fn', @alg_wso_compact_cec17); ...
    struct('name','BMJAYA', 'fn', @alg_bmjaya_budget_cec17); ...
    struct('name','DE'    , 'fn', @alg_de_cec17); ...
    struct('name','JAYA'  , 'fn', @alg_jaya_budget_cec17); ...
    struct('name','DOA'   , 'fn', @alg_doa_budget_cec17); ...
    struct('name','SBOA'  , 'fn', @alg_sboa_budget_cec17) ...
};

baseSeed = 20260602;

num_plot_points = 2001;
fes_show = unique(round(linspace(1, MaxFEs, num_plot_points)));

%% ============================= SETUP =================================
if ~exist(OUT_DIR, 'dir')
    mkdir(OUT_DIR);
end

addpath(genpath(CEC_DIR));
addpath(genpath(pwd));

fprintf('CEC2017 MEX found at:\n  %s\n\n', which('cec17_func'));

assert(exist('cec17_func', 'file') ~= 0, ...
    ['Cannot find cec17_func. Place the official CEC2017 MEX file and ', ...
     'input_data folder under CEC_DIR.']);

sanity_check_cec17();

%% ============================ STORAGE ================================
nDims = numel(Dims);
nAlgs = numel(Algorithms);
nFuncs = numel(FuncIDs);
nPts = numel(fes_show);
nCurveFuncs = numel(CurveFuncIDs);

MeanAll = nan(nDims, nFuncs, nAlgs);
StdAll = nan(nDims, nFuncs, nAlgs);
BestAll = nan(nDims, nFuncs, nAlgs);
MedianAll = nan(nDims, nFuncs, nAlgs);
WorstAll = nan(nDims, nFuncs, nAlgs);
RankAll = nan(nDims, nFuncs, nAlgs);

curve_error_sampled = nan(Runs, nCurveFuncs, nAlgs, nPts, 'single');

%% ============================= RUN ===================================
for d = 1:nDims
    D = Dims(d);
    [lb, ub] = cec17_bounds(D);

    for a = 1:nAlgs
        algo = Algorithms{a};

        out_csv = fullfile(OUT_DIR, ...
            sprintf('CEC2017_%s_%dD.csv', algo.name, D));

        fid = fopen(out_csv, 'w');
        assert(fid ~= -1, 'Cannot open CSV: %s', out_csv);

        fprintf(fid, ...
            'Function,Mean,Std,Best,Median,Worst,MaxFEs,PopSize,Runs,Dim\n');

        fprintf('\n=====================================================\n');
        fprintf('Algorithm: %s | Dim=%d | Pop=%d | MaxFEs=%d | Runs=%d\n', ...
            algo.name, D, N, MaxFEs, Runs);
        fprintf('Output: %s\n', out_csv);
        fprintf('=====================================================\n');

        for fpos = 1:nFuncs
            func_id = FuncIDs(fpos);

            [results, sampled_raw, fes_used] = run_one_function_cec17( ...
                algo.fn, func_id, D, N, MaxFEs, Runs, ...
                lb, ub, fes_show, baseSeed, a);

            assert(all(fes_used == MaxFEs), ...
                '%s F%d did not use exactly MaxFEs.', algo.name, func_id);

            mu = mean(results);
            sd = std(results);
            bst = min(results);
            med = median(results);
            wst = max(results);

            MeanAll(d, fpos, a) = mu;
            StdAll(d, fpos, a) = sd;
            BestAll(d, fpos, a) = bst;
            MedianAll(d, fpos, a) = med;
            WorstAll(d, fpos, a) = wst;

            fprintf('F%02d  mean=% .6e  std=% .6e  best=% .6e\n', ...
                func_id, mu, sd, bst);

            fprintf(fid, ...
                'F%d,%.15g,%.15g,%.15g,%.15g,%.15g,%d,%d,%d,%d\n', ...
                func_id, mu, sd, bst, med, wst, MaxFEs, N, Runs, D);

            curve_pos = find(CurveFuncIDs == func_id, 1);

            if ~isempty(curve_pos)
                fopt = 100 * func_id;

                curve_error_sampled(:, curve_pos, a, :) = single( ...
                    abs(double(sampled_raw) - fopt) + 1e-12);
            end
        end

        fclose(fid);
    end

    %% ======================== PER-FUNCTION RANK ======================
    for fpos = 1:nFuncs
        values = squeeze(MeanAll(d, fpos, :)).';

        RankAll(d, fpos, :) = reshape( ...
            average_rank_ascending(values), 1, 1, []);
    end

    %% ========================== EXPORT ===============================
    row_count = nFuncs * nAlgs;

    Function = cell(row_count, 1);
    Algorithm = cell(row_count, 1);

    Mean = nan(row_count, 1);
    Std = nan(row_count, 1);
    Best = nan(row_count, 1);
    Median = nan(row_count, 1);
    Worst = nan(row_count, 1);
    Rank = nan(row_count, 1);

    r = 0;

    for fpos = 1:nFuncs
        for a = 1:nAlgs
            r = r + 1;

            Function{r} = sprintf('F%d', FuncIDs(fpos));
            Algorithm{r} = Algorithms{a}.name;

            Mean(r) = MeanAll(d, fpos, a);
            Std(r) = StdAll(d, fpos, a);
            Best(r) = BestAll(d, fpos, a);
            Median(r) = MedianAll(d, fpos, a);
            Worst(r) = WorstAll(d, fpos, a);
            Rank(r) = RankAll(d, fpos, a);
        end
    end

    detailed_table = table( ...
        Function, Algorithm, Mean, Std, Best, Median, Worst, Rank);

    detailed_csv = fullfile(OUT_DIR, ...
        sprintf('CEC2017_ALL_Algorithms_%dD_with_Rank.csv', D));

    writetable(detailed_table, detailed_csv);

    MeanSci = cellfun(@(x) sprintf('%.3E', x), num2cell(Mean), ...
        'UniformOutput', false);
    StdSci = cellfun(@(x) sprintf('%.3E', x), num2cell(Std), ...
        'UniformOutput', false);
    BestSci = cellfun(@(x) sprintf('%.3E', x), num2cell(Best), ...
        'UniformOutput', false);
    MedianSci = cellfun(@(x) sprintf('%.3E', x), num2cell(Median), ...
        'UniformOutput', false);
    WorstSci = cellfun(@(x) sprintf('%.3E', x), num2cell(Worst), ...
        'UniformOutput', false);

    sci_table = table( ...
        Function, Algorithm, MeanSci, StdSci, BestSci, MedianSci, WorstSci, Rank, ...
        'VariableNames', { ...
            'Function','Algorithm','Mean','Std','Best','Median','Worst','Rank'});

    sci_csv = fullfile(OUT_DIR, ...
        sprintf('CEC2017_ALL_Algorithms_%dD_with_Rank_scientific.csv', D));

    writetable(sci_table, sci_csv);

    TotalRank = squeeze(sum(RankAll(d, :, :), 2));
    OverallRank = average_rank_ascending(TotalRank.');

    FinalAlgorithm = cell(nAlgs, 1);
    SumRank = nan(nAlgs, 1);
    FinalRank = nan(nAlgs, 1);

    for a = 1:nAlgs
        FinalAlgorithm{a} = Algorithms{a}.name;
        SumRank(a) = TotalRank(a);
        FinalRank(a) = OverallRank(a);
    end

    final_rank_table = table( ...
        FinalAlgorithm, SumRank, FinalRank, ...
        'VariableNames', {'Algorithm','SumRank','FinalRank'});

    final_rank_table = sortrows(final_rank_table, 'FinalRank');

    final_rank_csv = fullfile(OUT_DIR, ...
        sprintf('CEC2017_FinalRank_%dD.csv', D));

    writetable(final_rank_table, final_rank_csv);

    fprintf('\n=====================================================\n');
    fprintf('FINAL RANK TABLE | CEC2017 | Dim = %d\n', D);
    fprintf('Rule: smaller Mean(raw objective value) = better rank\n');
    fprintf('SumRank: sum of ranks over selected benchmark functions\n');
    fprintf('=====================================================\n');

    disp(final_rank_table);

    config = struct();
    config.alg_names = cellfun(@(s) s.name, Algorithms, 'UniformOutput', false);
    config.func_list = CurveFuncIDs;
    config.full_func_list = FuncIDs;
    config.fes_show = fes_show;
    config.D = D;
    config.N = N;
    config.MaxFEs = MaxFEs;
    config.Runs = Runs;
    config.curve_metric = '|f(x)-f*| + 1e-12';
    config.summary_metric = 'raw objective value f(x)';
    config.note = 'CEC2017 F2 excluded by default due to known instability.';

    result_mat = fullfile(OUT_DIR, ...
        sprintf('results_CEC2017_%dalg_D%d_N%d_FEs%d_Runs%d.mat', ...
            nAlgs, D, N, MaxFEs, Runs));

    save(result_mat, ...
        'config', 'MeanAll', 'StdAll', 'BestAll', ...
        'MedianAll', 'WorstAll', 'RankAll', ...
        'curve_error_sampled', 'final_rank_table', '-v7.3');

    fprintf('Detailed CSV:\n  %s\n', detailed_csv);
    fprintf('Scientific CSV:\n  %s\n', sci_csv);
    fprintf('Final Rank CSV:\n  %s\n', final_rank_csv);
    fprintf('MAT result:\n  %s\n', result_mat);
end

fprintf('\nAll done. Results saved to:\n  %s\n', OUT_DIR);
end

function sanity_check_cec17()
D = 50;
x = rand(D, 1) * 200 - 100;
f = cec17_func(x, 1);

fprintf('Sanity check: cec17_func(50x1) works. Example F1 value = %.6e\n\n', f);
end

function ranks = average_rank_ascending(values)
values = values(:).';

[sorted_values, order] = sort(values, 'ascend');

n = numel(values);
ranks = nan(1, n);

i = 1;

while i <= n
    j = i;

    tie_tol = 1e-12 * max(1, abs(sorted_values(i)));

    while j < n && ...
            abs(sorted_values(j + 1) - sorted_values(i)) <= tie_tol
        j = j + 1;
    end

    ranks(order(i:j)) = mean(i:j);
    i = j + 1;
end
end
