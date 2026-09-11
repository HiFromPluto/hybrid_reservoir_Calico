function R = report_narma_pack(varargin)
%REPORT_NARMA_PACK  Wringe pack 1 driver for the nine Narma10b dirs.
%
%   From the HybridDish folder (or anywhere on the MATLAB path):
%     addpath('matlab')
%     R = report_narma_pack
%
%   One driven seed only:
%     R = report_narma_pack('Only', "driven_seed222")
%
%   Uses frozen Narma10b voxels. Does not run BSim. Closed val-slice.
%   Stage 6 dirs are included as a second column when present.

    here = fileparts(mfilename('fullpath'));
    dish = fileparts(here);
    examples = fileparts(dish);
    narmaDir = fullfile(examples, 'BSimReservoirPlanNarma10b');
    stage6Dir = fullfile(examples, 'BSimReservoirPlanStage6');
    targetCsv = fullfile(narmaDir, 'narma10_target.csv');

    p = inputParser;
    addParameter(p, 'NarmaDir', narmaDir);
    addParameter(p, 'Stage6Dir', stage6Dir);
    addParameter(p, 'Target', targetCsv);
    addParameter(p, 'Only', "");
    addParameter(p, 'Horizons', [1, 5, 10]);
    parse(p, varargin{:});

    addpath(here);
    target = readtable(p.Results.Target, 'Delimiter', ';', 'FileType', 'text');
    y = target.y_next;
    u = target.u;
    B = baselines(y);

    arms = {'driven', 'brownian', 'silent'};
    seeds = [111, 222, 333];
    R.baselines = B;
    R.narma10b = struct();
    R.stage6 = struct();

    for a = 1:numel(arms)
        arm = arms{a};
        family = family_for(arm);
        for s = 1:numel(seeds)
            tag = sprintf('%s_seed%d', arm, seeds(s));
            if strlength(p.Results.Only) > 0 && ~strcmp(tag, char(p.Results.Only))
                continue
            end
            folder = fullfile(p.Results.NarmaDir, 'results', ['narma10b_' tag]);
            fprintf('\n=== Narma10b %s ===\n', tag);
            R.narma10b.(matlab.lang.makeValidName(tag)) = ...
                analyse_one(folder, family, y, u, p.Results.Horizons, arm);
        end
    end

    s6seeds = [101, 202, 303];
    for a = 1:numel(arms)
        arm = arms{a};
        family = family_for(arm);
        for s = 1:numel(s6seeds)
            tag = sprintf('%s_seed%d', arm, s6seeds(s));
            folder = fullfile(p.Results.Stage6Dir, 'results', ['stage6_' tag]);
            if ~isfile(fullfile(folder, 'voxels.csv'))
                continue
            end
            if strlength(p.Results.Only) > 0
                continue
            end
            fprintf('\n=== Stage 6, not this run %s ===\n', tag);
            R.stage6.(matlab.lang.makeValidName(tag)) = ...
                analyse_one(folder, family, y, u, p.Results.Horizons, arm);
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

function S = analyse_one(folder, family, y, u, horizons, arm)
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
    S.mc = memory_capacity(X, u);
    S.horizons = struct('h', {}, 'n_test', {}, 'test_nrmse', {}, 'test_r2', {}, 'lambda', {});
    for i = 1:numel(horizons)
        Hh = horizon_narma(X, y, horizons(i));
        S.horizons(i).h = Hh.h;
        S.horizons(i).n_test = Hh.n_test;
        S.horizons(i).test_nrmse = Hh.test_nrmse;
        S.horizons(i).test_r2 = Hh.test_r2;
        S.horizons(i).lambda = Hh.lambda;
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
