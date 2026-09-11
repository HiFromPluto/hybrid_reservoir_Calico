function R = report_mg_pack(varargin)
%REPORT_MG_PACK  Wringe pack 3 driver for BenchA Mackey-Glass dirs.
%
%   From the HybridDish folder:
%     addpath('matlab')
%     R = report_mg_pack
%
%   One driven seed (MATLAB reproduction target):
%     R = report_mg_pack('Only', "driven_seed101")
%
%   Frozen BenchA MG voxels. Does not run BSim. Closed val-slice.
%   Target is raw x_next, not affine u. Teacher-forced one-step, not
%   Jaeger autonomous MG. Python report_mg_pack.py writes WRINGE_MG_PACK.md.
%
%   Expected seed 101: driven NRMSE 0.3306, field-only 0.0229.

    here = fileparts(mfilename('fullpath'));
    dish = fileparts(here);
    examples = fileparts(dish);
    benchDir = fullfile(examples, 'BSimReservoirPlanBenchA');
    targetCsv = fullfile(benchDir, 'mg_target.csv');

    p = inputParser;
    addParameter(p, 'BenchDir', benchDir);
    addParameter(p, 'Target', targetCsv);
    addParameter(p, 'Only', "");
    addParameter(p, 'Horizons', [1, 5, 10]);
    parse(p, varargin{:});

    addpath(here);
    target = readtable(p.Results.Target, 'Delimiter', ';', 'FileType', 'text');
    y = target.x_next;
    B = baselines(y, 'Label', 'Mackey-Glass');

    arms = {'driven', 'brownian', 'silent'};
    seeds = [101, 202, 303];
    R.baselines = B;
    R.mg = struct();

    for a = 1:numel(arms)
        arm = arms{a};
        family = family_for(arm);
        for s = 1:numel(seeds)
            tag = sprintf('%s_seed%d', arm, seeds(s));
            if strlength(p.Results.Only) > 0 && ~strcmp(tag, char(p.Results.Only))
                continue
            end
            folder = fullfile(p.Results.BenchDir, 'results', ['mg_' tag]);
            fprintf('\n=== Mackey-Glass %s ===\n', tag);
            R.mg.(matlab.lang.makeValidName(tag)) = ...
                analyse_one(folder, family, y, p.Results.Horizons, arm);
        end
    end
end

function family = family_for(arm)
    if strcmp(arm, 'brownian')
        family = 'brownian';
    else
        family = 'biology';
    end
end

function S = analyse_one(folder, family, y, horizons, arm)
    T = load_voxels(folder);
    [X, names] = window_features(T, family);
    fprintf('  family=%s  n_features=%d\n', family, numel(names));
    h1 = ridge_closed(X, y);
    fprintf('  h=1 NRMSE=%.4f  R^2=%.4f  lambda=%.4g\n', ...
        h1.test_nrmse, h1.test_r2, h1.lambda);
    S.folder = folder;
    S.family = family;
    S.n_features = numel(names);
    S.h1 = h1;
    S.horizons = struct('h', {}, 'n_test', {}, 'test_nrmse', {}, 'test_r2', {}, 'lambda', {});
    if strcmp(arm, 'driven')
        for i = 1:numel(horizons)
            Hh = horizon_narma(X, y, horizons(i));
            S.horizons(i).h = Hh.h;
            S.horizons(i).n_test = Hh.n_test;
            S.horizons(i).test_nrmse = Hh.test_nrmse;
            S.horizons(i).test_r2 = Hh.test_r2;
            S.horizons(i).lambda = Hh.lambda;
        end
    end
    S.field = [];
    if strcmp(arm, 'driven')
        [Xf, ~] = window_features(T, 'field');
        Sf = ridge_closed(Xf, y);
        fprintf('  field-only h=1 NRMSE=%.4f\n', Sf.test_nrmse);
        S.field = Sf;
        S.field_horizons = struct('h', {}, 'n_test', {}, 'test_nrmse', {}, 'test_r2', {}, 'lambda', {});
        for i = 1:numel(horizons)
            Fh = horizon_narma(Xf, y, horizons(i));
            S.field_horizons(i).h = Fh.h;
            S.field_horizons(i).n_test = Fh.n_test;
            S.field_horizons(i).test_nrmse = Fh.test_nrmse;
            S.field_horizons(i).test_r2 = Fh.test_r2;
            S.field_horizons(i).lambda = Fh.lambda;
        end
    end
end
