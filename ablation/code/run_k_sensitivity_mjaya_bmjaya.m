function run_k_sensitivity_mjaya_bmjaya()
clc; clear; close all;

% Sensitivity test for the number of subpopulations K.
% Algorithms:
%   MJAYA  = multipopulation JAYA + elite sharing
%   BMJAYA = BMR repair + multipopulation evolution + elite sharing

SCRIPT_DIR = fileparts(mfilename('fullpath'));
ROOT_DIR = fileparts(SCRIPT_DIR);
OUT_DIR = fullfile(ROOT_DIR, 'results_k_sensitivity');

Dims = 20;           % change to [10, 20] if both dimensions are needed
N = 50;
MaxFEs = 1e5;
Runs = 30;
FuncIDs = 1:12;
KList = 1:5;

Algorithms = { ...
    struct('name', 'MJAYA',  'fn', @alg_mjaya_fixedN_budget); ...
    struct('name', 'BMJAYA', 'fn', @alg_bmjaya_fixedN_budget) ...
};

if ~exist(OUT_DIR, 'dir')
    mkdir(OUT_DIR);
end

addpath(genpath(SCRIPT_DIR));
fprintf("CEC2022 MEX found at:\n  %s\n\n", which('cec22_test_func'));
sanity_check_cec22();

for d = 1:numel(Dims)
    D = Dims(d);
    [lb, ub] = cec_bounds(D);

    SummaryRows = {};
    DetailRows = {};

    for a = 1:numel(Algorithms)
        algo = Algorithms{a};

        for K = KList
            out_csv = fullfile(OUT_DIR, sprintf('CEC2022_%s_K%d_%dD.csv', algo.name, K, D));
            fid = fopen(out_csv, 'w');
            assert(fid ~= -1, 'Cannot open output CSV: %s', out_csv);
            fprintf(fid, "Function,Mean,Std,Best,Median,Worst,MaxFEs,PopSize,Runs,Dim,K\n");

            fprintf("\n=====================================================\n");
            fprintf("Algorithm: %s | K=%d | Dim=%d | Pop=%d | MaxFEs=%d | Runs=%d\n", ...
                algo.name, K, D, N, MaxFEs, Runs);
            fprintf("Output: %s\n", out_csv);
            fprintf("=====================================================\n");

            for func_id = FuncIDs
                results = zeros(Runs, 1);

                for r = 1:Runs
                    rng(1 + r, 'twister');
                    results(r) = algo.fn(func_id, D, N, MaxFEs, lb, ub, K);
                end

                mu = mean(results);
                sd = std(results);
                bst = min(results);
                med = median(results);
                wst = max(results);

                fprintf("F%02d  mean=% .6e  std=% .6e  best=% .6e\n", func_id, mu, sd, bst);
                fprintf(fid, "F%d,%.15g,%.15g,%.15g,%.15g,%.15g,%d,%d,%d,%d,%d\n", ...
                    func_id, mu, sd, bst, med, wst, MaxFEs, N, Runs, D, K);

                DetailRows(end+1, :) = {D, algo.name, K, func_id, mu, sd, bst, med, wst}; %#ok<AGROW>
            end

            fclose(fid);
        end
    end

    Detail = cell2table(DetailRows, 'VariableNames', ...
        {'Dim','Algorithm','K','Function','Mean','Std','Best','Median','Worst'});

    Detail.Rank = zeros(height(Detail), 1);
    for func_id = FuncIDs
        idx = Detail.Function == func_id;
        values = Detail.Mean(idx);
        Detail.Rank(idx) = average_rank_ascending(values);
    end

    groups = unique(Detail(:, {'Algorithm', 'K'}), 'rows');
    for i = 1:height(groups)
        idx = strcmp(Detail.Algorithm, groups.Algorithm{i}) & Detail.K == groups.K(i);
        total_rank = sum(Detail.Rank(idx));
        avg_rank = mean(Detail.Rank(idx));
        mean_value = mean(Detail.Mean(idx));
        SummaryRows(end+1, :) = {D, groups.Algorithm{i}, groups.K(i), total_rank, avg_rank, mean_value}; %#ok<AGROW>
    end

    Summary = cell2table(SummaryRows, 'VariableNames', ...
        {'Dim','Algorithm','K','TotalRank','AverageRank','MeanOfMeans'});
    Summary = sortrows(Summary, {'TotalRank', 'AverageRank'});

    detail_csv = fullfile(OUT_DIR, sprintf('CEC2022_MJAYA_BMJAYA_K1to5_D%d_detail.csv', D));
    summary_csv = fullfile(OUT_DIR, sprintf('CEC2022_MJAYA_BMJAYA_K1to5_D%d_summary.csv', D));
    writetable(Detail, detail_csv);
    writetable(Summary, summary_csv);

    fprintf("\n=====================================================\n");
    fprintf("K sensitivity summary | D=%d\n", D);
    fprintf("=====================================================\n");
    disp(Summary);
    fprintf("Detail CSV:  %s\n", detail_csv);
    fprintf("Summary CSV: %s\n", summary_csv);
end

fprintf("\nAll done. CSV files saved to: %s\n", OUT_DIR);
end

function ranks = average_rank_ascending(values)
[sorted_values, order] = sort(values(:), 'ascend');
ranks = zeros(size(values(:)));
i = 1;
n = numel(values);

while i <= n
    j = i;
    while j < n && sorted_values(j + 1) == sorted_values(i)
        j = j + 1;
    end

    avg_rank = (i + j) / 2;
    ranks(order(i:j)) = avg_rank;
    i = j + 1;
end
end

function sanity_check_cec22()
D = 20;
x = rand(D, 1) * 200 - 100;
f = cec22_test_func(x, 7);
fprintf("Sanity check: cec22_test_func(Dx1) works. Example F7 value = %.6e\n\n", f);
end
