function run_cec2022_all()
clc; clear; close all;

% =========================
% User config
% =========================
CEC_DIR   = fileparts(mfilename('fullpath'));  % directory containing cec22_test_func.mexw64 and input_data/
OUT_DIR   = fullfile(pwd, 'results');   % output folder

Dims      = [10, 20];     % run both 10D and 20D
N         = 50;           % population size
MaxFEs    = 1e5;          % strict budget
Runs      = 30;           % independent runs
FuncIDs   = 1:12;         % CEC2022 F1..F12

Algorithms = { ...
    struct('name','PSO'   , 'fn', @alg_pso_weak); ...
    struct('name','GWO'   , 'fn', @alg_gwo); ...
    struct('name','WSO'   , 'fn', @alg_wso_compact); ...
    struct('name','BMJAYA', 'fn', @alg_bmjaya_wrapper); ...
    struct('name','DE'    , 'fn', @alg_de); ...
    struct('name','JAYA'  , 'fn', @alg_jaya_budget); ...
    struct('name','SBOA'  , 'fn', @alg_sboa_wrapper); ...
    struct('name','DOA'   , 'fn', @alg_doa_wrapper) ...
};

% =========================
% Setup
% =========================
if ~exist(OUT_DIR, 'dir')
    mkdir(OUT_DIR);
end

addpath(genpath(CEC_DIR));      % make mex + input_data visible
addpath(genpath(pwd));          % make this project visible

fprintf("CEC2022 MEX found at:\n  %s\n\n", which('cec22_test_func'));

% quick sanity check
sanity_check_cec22();

% =========================
% Rank storage
% =========================
nDims  = numel(Dims);
nAlgs  = numel(Algorithms);
nFuncs = numel(FuncIDs);

% Dimensions:
%   dim_index x function_position x algorithm_index
MeanAll   = nan(nDims, nFuncs, nAlgs);
StdAll    = nan(nDims, nFuncs, nAlgs);
BestAll   = nan(nDims, nFuncs, nAlgs);
MedianAll = nan(nDims, nFuncs, nAlgs);
WorstAll  = nan(nDims, nFuncs, nAlgs);
RankAll   = nan(nDims, nFuncs, nAlgs);

% Runtime storage
TimeTotalAll = nan(nDims, nFuncs, nAlgs);   % total time for one function over Runs
TimeMeanAll  = nan(nDims, nFuncs, nAlgs);   % average time per run
AlgTimeAll   = nan(nDims, nAlgs);           % total time for one algorithm over all functions

% =========================
% Run all
% =========================
total_tic = tic;

for d = 1:nDims
    D = Dims(d);
    [lb, ub] = cec_bounds(D);

    fprintf("\n\n#####################################################\n");
    fprintf("START DIMENSION: D = %d\n", D);
    fprintf("#####################################################\n");

    dim_tic = tic;

    for a = 1:nAlgs
        algo = Algorithms{a};

        out_csv = fullfile(OUT_DIR, sprintf('CEC2022_%s_%dD.csv', algo.name, D));
        fid = fopen(out_csv, 'w');

        assert(fid ~= -1, 'Cannot open output CSV: %s', out_csv);

        fprintf(fid, ...
            "Function,Mean,Std,Best,Median,Worst,TimeTotal_s,TimeMean_s,MaxFEs,PopSize,Runs,Dim\n");

        fprintf("\n=====================================================\n");
        fprintf("Algorithm: %s | Dim=%d | Pop=%d | MaxFEs=%d | Runs=%d\n", ...
            algo.name, D, N, MaxFEs, Runs);
        fprintf("Output: %s\n", out_csv);
        fprintf("=====================================================\n");

        alg_tic = tic;

        for fpos = 1:nFuncs
            func_id = FuncIDs(fpos);

            func_tic = tic;

            results = run_one_function( ...
                algo.fn, func_id, D, N, MaxFEs, Runs, lb, ub);

            func_time_total = toc(func_tic);
            func_time_mean  = func_time_total / Runs;

            mu  = mean(results);
            sd  = std(results);
            bst = min(results);
            med = median(results);
            wst = max(results);

            % Store raw objective statistics f(x).
            MeanAll(d, fpos, a)   = mu;
            StdAll(d, fpos, a)    = sd;
            BestAll(d, fpos, a)   = bst;
            MedianAll(d, fpos, a) = med;
            WorstAll(d, fpos, a)  = wst;

            % Store runtime.
            TimeTotalAll(d, fpos, a) = func_time_total;
            TimeMeanAll(d, fpos, a)  = func_time_mean;

            fprintf("F%02d  mean=% .6e  std=% .6e  best=% .6e  time=%.2fs  avg=%.2fs/run\n", ...
                func_id, mu, sd, bst, func_time_total, func_time_mean);

            fprintf(fid, ...
                "F%d,%.15g,%.15g,%.15g,%.15g,%.15g,%.6f,%.6f,%d,%d,%d,%d\n", ...
                func_id, mu, sd, bst, med, wst, ...
                func_time_total, func_time_mean, MaxFEs, N, Runs, D);
        end

        alg_time_total = toc(alg_tic);
        AlgTimeAll(d, a) = alg_time_total;

        fprintf("\nAlgorithm %s | Dim=%d finished.\n", algo.name, D);
        fprintf("Total algorithm time = %.2f s = %.2f min = %.2f h\n", ...
            alg_time_total, alg_time_total / 60, alg_time_total / 3600);

        fclose(fid);
    end

    dim_time_total = toc(dim_tic);

    fprintf("\n#####################################################\n");
    fprintf("DIMENSION D=%d finished. Total time = %.2f s = %.2f min = %.2f h\n", ...
        D, dim_time_total, dim_time_total / 60, dim_time_total / 3600);
    fprintf("#####################################################\n");

    % =========================
    % Calculate per-function rank
    % =========================
    for fpos = 1:nFuncs
        values = squeeze(MeanAll(d, fpos, :)).';
        RankAll(d, fpos, :) = reshape( ...
            average_rank_ascending(values), 1, 1, []);
    end

    % =========================
    % Export combined detailed table
    % =========================
    row_count = nFuncs * nAlgs;

    Function = cell(row_count, 1);
    Algorithm = cell(row_count, 1);
    Mean = nan(row_count, 1);
    Std = nan(row_count, 1);
    Best = nan(row_count, 1);
    Median = nan(row_count, 1);
    Worst = nan(row_count, 1);
    TimeTotal_s = nan(row_count, 1);
    TimeMean_s = nan(row_count, 1);
    Rank = nan(row_count, 1);

    r = 0;
    for fpos = 1:nFuncs
        for a = 1:nAlgs
            r = r + 1;

            Function{r}  = sprintf('F%d', FuncIDs(fpos));
            Algorithm{r} = Algorithms{a}.name;

            Mean(r)   = MeanAll(d, fpos, a);
            Std(r)    = StdAll(d, fpos, a);
            Best(r)   = BestAll(d, fpos, a);
            Median(r) = MedianAll(d, fpos, a);
            Worst(r)  = WorstAll(d, fpos, a);

            TimeTotal_s(r) = TimeTotalAll(d, fpos, a);
            TimeMean_s(r)  = TimeMeanAll(d, fpos, a);

            Rank(r) = RankAll(d, fpos, a);
        end
    end

    detailed_rank_table = table( ...
        Function, Algorithm, Mean, Std, Best, Median, Worst, ...
        TimeTotal_s, TimeMean_s, Rank);

    detailed_rank_csv = fullfile(OUT_DIR, ...
        sprintf('CEC2022_ALL_Algorithms_%dD_with_Rank.csv', D));

    writetable(detailed_rank_table, detailed_rank_csv);

    % =========================
    % Export scientific-notation detailed table
    % =========================
    MeanSci   = cell(row_count, 1);
    StdSci    = cell(row_count, 1);
    BestSci   = cell(row_count, 1);
    MedianSci = cell(row_count, 1);
    WorstSci  = cell(row_count, 1);

    TimeTotalSci = cell(row_count, 1);
    TimeMeanSci  = cell(row_count, 1);

    for r = 1:row_count
        MeanSci{r}   = sprintf('%.3E', Mean(r));
        StdSci{r}    = sprintf('%.3E', Std(r));
        BestSci{r}   = sprintf('%.3E', Best(r));
        MedianSci{r} = sprintf('%.3E', Median(r));
        WorstSci{r}  = sprintf('%.3E', Worst(r));

        TimeTotalSci{r} = sprintf('%.3E', TimeTotal_s(r));
        TimeMeanSci{r}  = sprintf('%.3E', TimeMean_s(r));
    end

    detailed_rank_sci_table = table( ...
        Function, Algorithm, MeanSci, StdSci, BestSci, ...
        MedianSci, WorstSci, TimeTotalSci, TimeMeanSci, Rank, ...
        'VariableNames', { ...
            'Function','Algorithm','Mean','Std','Best', ...
            'Median','Worst','TimeTotal_s','TimeMean_s','Rank'});

    detailed_rank_sci_csv = fullfile(OUT_DIR, ...
        sprintf('CEC2022_ALL_Algorithms_%dD_with_Rank_scientific.csv', D));

    writetable(detailed_rank_sci_table, detailed_rank_sci_csv);

    % =========================
    % Calculate final overall rank
    % =========================
    % Sum F1-F12 ranks for every algorithm.
    TotalRank = squeeze(sum(RankAll(d, :, :), 2));
    OverallRank = average_rank_ascending(TotalRank.');

    Algorithm = cell(nAlgs, 1);
    SumRank = nan(nAlgs, 1);
    FinalRank = nan(nAlgs, 1);
    AlgTimeTotal_s = nan(nAlgs, 1);
    AlgTimeMeanPerFunc_s = nan(nAlgs, 1);

    for a = 1:nAlgs
        Algorithm{a} = Algorithms{a}.name;
        SumRank(a) = TotalRank(a);
        FinalRank(a) = OverallRank(a);

        AlgTimeTotal_s(a) = AlgTimeAll(d, a);
        AlgTimeMeanPerFunc_s(a) = AlgTimeAll(d, a) / nFuncs;
    end

    final_rank_table = table( ...
        Algorithm, SumRank, FinalRank, AlgTimeTotal_s, AlgTimeMeanPerFunc_s);

    final_rank_table = sortrows(final_rank_table, 'FinalRank');

    final_rank_csv = fullfile(OUT_DIR, ...
        sprintf('CEC2022_FinalRank_%dD.csv', D));

    writetable(final_rank_table, final_rank_csv);

    % =========================
    % Print rank results
    % =========================
    fprintf("\n=====================================================\n");
    fprintf("FINAL RANK TABLE | Dim = %d\n", D);
    fprintf("Rule: smaller Mean(raw objective value) = better rank\n");
    fprintf("SumRank: sum of F1-F12 ranks\n");
    fprintf("Runtime: measured by tic/toc in seconds\n");
    fprintf("=====================================================\n");

    disp(final_rank_table);

    fprintf("Detailed numeric Rank CSV:\n  %s\n", detailed_rank_csv);
    fprintf("Detailed scientific Rank CSV:\n  %s\n", detailed_rank_sci_csv);
    fprintf("Final Rank CSV:\n  %s\n", final_rank_csv);
end

total_time = toc(total_tic);

fprintf("\nAll done. CSV files saved to: %s\n", OUT_DIR);
fprintf("Total running time = %.2f s = %.2f min = %.2f h\n", ...
    total_time, total_time / 60, total_time / 3600);

end

% -------------------------
function sanity_check_cec22()
D = 20;
x = rand(D,1) * 200 - 100;  % D脳1 required by mex
f = cec22_test_func(x, 7);

fprintf("Sanity check: cec22_test_func(D脳1) works. Example F7 value = %.6e\n\n", f);
end

% -------------------------
function ranks = average_rank_ascending(values)
% AVERAGE_RANK_ASCENDING
% Smaller values receive better ranks.
% Tied values receive their average rank.
% No Statistics Toolbox is required.

values = values(:).';

[sorted_values, order] = sort(values, 'ascend');

n = numel(values);
ranks = nan(1, n);

i = 1;

while i <= n
    j = i;

    % Numerical tolerance for tie detection.
    tie_tol = 1e-12 * max(1, abs(sorted_values(i)));

    while j < n && ...
            abs(sorted_values(j + 1) - sorted_values(i)) <= tie_tol
        j = j + 1;
    end

    ranks(order(i:j)) = mean(i:j);

    i = j + 1;
end
end
