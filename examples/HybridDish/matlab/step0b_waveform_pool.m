function T = step0b_waveform_pool(varargin)
%STEP0B_WAVEFORM_POOL  Coarser pooling of frozen 20x10 waveform maps.
%   classify_waveform rules (val AUC, not NRMSE). Does not run BSim.

    here = fileparts(mfilename('fullpath'));
    addpath(here);
    examples = fileparts(fileparts(here));
    waveDir = fullfile(examples, 'BSimReservoirPlanWaveform');
    outDir = fullfile(waveDir, 'results');
    labelsCsv = fullfile(waveDir, 'waveform_labels.csv');

    p = inputParser;
    addParameter(p, 'WaveDir', waveDir);
    addParameter(p, 'OutDir', outDir);
    addParameter(p, 'Labels', labelsCsv);
    parse(p, varargin{:});
    waveDir = p.Results.WaveDir;
    outDir = p.Results.OutDir;

    L = readtable(p.Results.Labels, 'Delimiter', ';', 'FileType', 'text');
    y = L.y(:);

    pools = {
        'P0', 1, 1, 20, 10
        'P1', 2, 2, 10, 5
        'P2', 4, 5, 5, 2
        'P3', 20, 10, 1, 1
        };
    seeds = [101, 202, 303];
    skipped = {};
    rows = {};

    for s = 1:numel(seeds)
        seed = seeds(s);
        folder = fullfile(waveDir, 'results', sprintf('waveform_driven_seed%d', seed));
        csvPath = fullfile(folder, 'voxels.csv');
        if ~isfile(csvPath)
            fprintf('SKIP missing %s\n', csvPath);
            skipped{end+1} = csvPath; %#ok<AGROW>
            continue
        end
        fprintf('\n=== waveform driven seed %d ===\n', seed);
        Tab = load_voxels(folder);
        X408 = window_features(Tab, 'biology');
        XR = window_features(Tab, 'receiver');
        XL = window_features(Tab, 'lum');
        XF = window_features(Tab, 'field');
        if size(X408, 2) ~= 408 || size(XR, 2) ~= 200 || size(XL, 2) ~= 200 || size(XF, 2) ~= 200
            error('step0b:cols', 'Unexpected map widths seed %d', seed);
        end
        deaths = X408(:, 401:408);
        rows = [rows; fit_driven_pools(seed, XR, XL, XF, deaths, y, pools)]; %#ok<AGROW>
    end

    for s = 1:numel(seeds)
        seed = seeds(s);
        folder = fullfile(waveDir, 'results', sprintf('waveform_brownian_seed%d', seed));
        csvPath = fullfile(folder, 'voxels.csv');
        if ~isfile(csvPath)
            fprintf('SKIP missing %s\n', csvPath);
            skipped{end+1} = csvPath; %#ok<AGROW>
            continue
        end
        fprintf('\n=== waveform brownian seed %d ===\n', seed);
        Tab = load_voxels(folder);
        XD = window_features(Tab, 'brownian');
        if size(XD, 2) ~= 200
            error('step0b:cols', 'Brownian Den width %d', size(XD, 2));
        end
        for ip = 1:size(pools, 1)
            tag = pools{ip, 1};
            bx = pools{ip, 2};
            by = pools{ip, 3};
            nx = pools{ip, 4};
            ny = pools{ip, 5};
            Dp = pool_state_map(XD, bx, by);
            nFeat = size(Dp, 2);
            if nFeat ~= nx * ny
                error('step0b:nfeat', 'Fden %s n_features=%d', tag, nFeat);
            end
            fit = classify_waveform_features(Dp, y, 'Quiet', true);
            fprintf('  %s Fden n=%d lambda=%.4g AUC=%.4f acc=%.4f\n', ...
                tag, nFeat, fit.lambda, fit.macro_auc, fit.accuracy);
            rows(end+1, :) = {'waveform', seed, 'brownian', 'Fden', tag, nFeat, ...
                fit.lambda, fit.macro_auc, fit.accuracy}; %#ok<AGROW>
        end
    end

    silentFolder = fullfile(waveDir, 'results', 'waveform_silent_seed101');
    if isfile(fullfile(silentFolder, 'voxels.csv'))
        fprintf('\n=== waveform silent seed 101 (F408 P0 sanity only) ===\n');
        Tab = load_voxels(silentFolder);
        X408 = window_features(Tab, 'biology');
        fit = classify_waveform_features(X408, y, 'Quiet', true);
        fprintf('  P0 F408 n=408 lambda=%.4g AUC=%.4f acc=%.4f\n', ...
            fit.lambda, fit.macro_auc, fit.accuracy);
        rows(end+1, :) = {'waveform', 101, 'silent', 'F408', 'P0', 408, ...
            fit.lambda, fit.macro_auc, fit.accuracy};
    else
        fprintf('SKIP missing silent seed 101 voxels (optional sanity)\n');
        skipped{end+1} = fullfile(silentFolder, 'voxels.csv');
    end

    T = cell2table(rows, 'VariableNames', {
        'task', 'seed', 'arm', 'family', 'pool', 'n_features', ...
        'lambda', 'test_macro_auc', 'test_acc'});
    if ~exist(outDir, 'dir')
        mkdir(outDir);
    end
    outCsv = fullfile(outDir, 'step0b_waveform_pool.csv');
    writetable(T, outCsv);
    fprintf('\nWrote %s (%d rows)\n', outCsv, height(T));
    if ~isempty(skipped)
        fprintf('Skipped %d missing path(s)\n', numel(skipped));
    end
    check_waveform_sanity(T);
end

function rows = fit_driven_pools(seed, XR, XL, XF, deaths, y, pools)
    rows = {};
    for ip = 1:size(pools, 1)
        tag = pools{ip, 1};
        bx = pools{ip, 2};
        by = pools{ip, 3};
        nx = pools{ip, 4};
        ny = pools{ip, 5};
        Rp = pool_state_map(XR, bx, by);
        Lp = pool_state_map(XL, bx, by);
        Fp = pool_state_map(XF, bx, by);
        cmap = size(Rp, 2);
        if cmap ~= nx * ny
            error('step0b:poolcount', '%s expected %d cols, got %d', tag, nx*ny, cmap);
        end
        families = {
            'F408', [Rp, Lp, deaths], 2*cmap + 8
            'FRL',  [Rp, Lp],         2*cmap
            'FR',   Rp,               cmap
            'FL',   Lp,               cmap
            'Ffield', Fp,             cmap
            };
        for f = 1:size(families, 1)
            fam = families{f, 1};
            X = families{f, 2};
            nFeat = families{f, 3};
            if size(X, 2) ~= nFeat
                error('step0b:nfeat', '%s %s n_features=%d expected %d', ...
                    tag, fam, size(X, 2), nFeat);
            end
            fit = classify_waveform_features(X, y, 'Quiet', true);
            fprintf('  %s %s n=%d lambda=%.4g AUC=%.4f acc=%.4f\n', ...
                tag, fam, nFeat, fit.lambda, fit.macro_auc, fit.accuracy);
            rows(end+1, :) = {'waveform', seed, 'driven', fam, tag, nFeat, ...
                fit.lambda, fit.macro_auc, fit.accuracy}; %#ok<AGROW>
        end
    end
end

function check_waveform_sanity(T)
    expAuc = [0.7969, 0.8287, 0.8411];
    expAcc = [0.76, 0.76, 0.68];
    seeds = [101, 202, 303];
    for i = 1:3
        got = pick(T, seeds(i), 'driven', 'F408', 'P0', 'test_macro_auc');
        if isempty(got)
            error('step0b:sanity', 'STOP missing driven F408 P0 seed %d', seeds(i));
        end
        if abs(got - expAuc(i)) > 1e-3
            error('step0b:sanity', ...
                'STOP F408 P0 seed %d AUC=%.6f expected %.4f', seeds(i), got, expAuc(i));
        end
        acc = pick(T, seeds(i), 'driven', 'F408', 'P0', 'test_acc');
        if abs(acc - expAcc(i)) > 1e-3
            error('step0b:sanity', ...
                'STOP F408 P0 seed %d acc=%.6f expected %.2f', seeds(i), acc, expAcc(i));
        end
        fld = pick(T, seeds(i), 'driven', 'Ffield', 'P0', 'test_macro_auc');
        if abs(fld - 0.8603) > 1e-3
            error('step0b:sanity', ...
                'STOP Ffield P0 seed %d AUC=%.6f expected 0.8603', seeds(i), fld);
        end
        facc = pick(T, seeds(i), 'driven', 'Ffield', 'P0', 'test_acc');
        if abs(facc - 0.86) > 1e-3
            error('step0b:sanity', ...
                'STOP Ffield P0 seed %d acc=%.6f expected 0.86', seeds(i), facc);
        end
    end
    denAuc = [];
    denAcc = [];
    for i = 1:3
        a = pick(T, seeds(i), 'brownian', 'Fden', 'P0', 'test_macro_auc');
        b = pick(T, seeds(i), 'brownian', 'Fden', 'P0', 'test_acc');
        if isempty(a)
            error('step0b:sanity', 'STOP missing brownian Fden P0 seed %d', seeds(i));
        end
        denAuc(end+1) = a; %#ok<AGROW>
        denAcc(end+1) = b; %#ok<AGROW>
    end
    if abs(mean(denAuc) - 0.5867) > 1e-3
        error('step0b:sanity', 'STOP Fden P0 mean AUC=%.6f expected ~0.5867', mean(denAuc));
    end
    if abs(mean(denAcc) - 0.31) > 2e-2
        error('step0b:sanity', 'STOP Fden P0 mean acc=%.6f expected ~0.31', mean(denAcc));
    end
    silent = pick(T, 101, 'silent', 'F408', 'P0', 'test_macro_auc');
    if ~isempty(silent) && abs(silent - 0.5) > 2e-2
        error('step0b:sanity', 'STOP silent F408 P0 AUC=%.6f expected ~0.5', silent);
    end
    fprintf('SANITY PASS  F408 AUC 0.7969/0.8287/0.8411  field 0.8603  Den ~0.5867\n');
end

function v = pick(T, seed, arm, family, pool, col)
    hit = T.seed == seed & strcmp(T.arm, arm) & strcmp(T.family, family) & strcmp(T.pool, pool);
    if sum(hit) == 0
        v = [];
        return
    end
    if sum(hit) ~= 1
        error('step0b:pick', 'expected 1 row for %s %s %s %d', arm, family, pool, seed);
    end
    v = T.(col)(hit);
end
