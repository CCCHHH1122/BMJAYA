function run_cec2022_all()
clc; clear; close all;

% =========================
% User config
% =========================
SCRIPT_DIR = fileparts(mfilename('fullpath'));
ROOT_DIR = fileparts(SCRIPT_DIR);
CEC_DIR   = SCRIPT_DIR;                       % contains cec22_test_func.mexw64 and input_data/
OUT_DIR   = fullfile(ROOT_DIR, 'results_ablation_20D');  % output folder

Dims      = 20;           % Table 2 ablation setting
N         = 50;           % population size
MaxFEs    = 1e5;          % strict budget
Runs      = 30;           % independent runs
FuncIDs   = 1:12;         % CEC2022 F1..F12

% Ablation variants:
% JAYA   = original JAYA
% BJAYA  = JAYA + BMR repair
% MJAYA  = multipopulation JAYA + elite sharing
% BMJAYA = BJAYA + multipopulation evolution + elite sharing
Algorithms = { ...
    struct('name','JAYA',   'fn', @alg_jaya_budget); ...
    struct('name','BJAYA',  'fn', @alg_bjaya_budget); ...
    struct('name','MJAYA',  'fn', @alg_mjaya_wrapper); ...
    struct('name','BMJAYA', 'fn', @alg_bmjaya_wrapper) ...
};

% =========================
% Setup
% =========================
if ~exist(OUT_DIR, 'dir')
    mkdir(OUT_DIR);
end

addpath(genpath(CEC_DIR));
fprintf("CEC2022 MEX found at:\n  %s\n\n", which('cec22_test_func'));
sanity_check_cec22();

% =========================
% Run all
% =========================
for d = 1:numel(Dims)
    D = Dims(d);
    [lb, ub] = cec_bounds(D);

    for a = 1:numel(Algorithms)
        algo = Algorithms{a};

        out_csv = fullfile(OUT_DIR, sprintf('CEC2022_%s_%dD.csv', algo.name, D));
        fid = fopen(out_csv, 'w');
        assert(fid ~= -1, 'Cannot open output CSV: %s', out_csv);

        fprintf(fid, "Function,Mean,Std,Best,Median,Worst,MaxFEs,PopSize,Runs,Dim\n");

        fprintf("\n=====================================================\n");
        fprintf("Algorithm: %s | Dim=%d | Pop=%d | MaxFEs=%d | Runs=%d\n", ...
            algo.name, D, N, MaxFEs, Runs);
        fprintf("Output: %s\n", out_csv);
        fprintf("=====================================================\n");

        for func_id = FuncIDs
            results = run_one_function(algo.fn, func_id, D, N, MaxFEs, Runs, lb, ub);

            mu  = mean(results);
            sd  = std(results);
            bst = min(results);
            med = median(results);
            wst = max(results);

            fprintf("F%02d  mean=% .6e  std=% .6e  best=% .6e\n", func_id, mu, sd, bst);

            fprintf(fid, "F%d,%.15g,%.15g,%.15g,%.15g,%.15g,%d,%d,%d,%d\n", ...
                func_id, mu, sd, bst, med, wst, MaxFEs, N, Runs, D);
        end

        fclose(fid);
    end
end

fprintf("\nAll done. CSV files saved to: %s\n", OUT_DIR);
end

function sanity_check_cec22()
D = 20;
x = rand(D, 1) * 200 - 100;
f = cec22_test_func(x, 7);
fprintf("Sanity check: cec22_test_func(Dx1) works. Example F7 value = %.6e\n\n", f);
end
