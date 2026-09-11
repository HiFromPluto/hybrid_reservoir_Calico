function T = step0b_mg_pool(varargin)
%STEP0B_MG_POOL  Coarser pooling of frozen 20x10 Mackey-Glass maps.
%   ridge_closed unchanged. Teacher-forced h=1 only. Does not run BSim.

    here = fileparts(mfilename('fullpath'));
    addpath(here);
    examples = fileparts(fileparts(here));
    benchDir = fullfile(examples, 'BSimReservoirPlanBenchA');
    outDir = fullfile(benchDir, 'results');
    targetCsv = fullfile(benchDir, 'mg_target.csv');

    p = inputParser;
    addParameter(p, 'BenchDir', benchDir);
    addParameter(p, 'OutDir', outDir);
    addParameter(p, 'Target', targetCsv);
    parse(p, varargin{:});
    benchDir = p.Results.BenchDir;
    outDir = p.Results.OutDir;

    target = readtable(p.Results.Target, 'Delimiter', ';', 'FileType', 'text');
    y = target.x_next(:);
    B = baselines(y, 'Label', 'Mackey-Glass');
    fprintf('  persistence=%.4f  intercept=%.4f\n', ...
        B.persistence_nrmse, B.intercept_nrmse);

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
        folder = fullfile(benchDir, 'results', sprintf('mg_driven_seed%d', seed));
        csvPath = fullfile(folder, 'voxels.csv');
        if ~isfile(csvPath)
            fprintf('SKIP missing %s\n', csvPath);
            skipped{end+1} = csvPath; %#ok<AGROW>
            continue
        end
        fprintf('\n=== MG driven seed %d ===\n', seed);
        Tab = load_voxels(folder);
        X408 = window_features(Tab, 'biology');
        XR = window_features(Tab, 'receiver');
        XL = window_features(Tab, 'lum');
        XF = window_features(Tab, 'field');
        if size(X408, 2) ~= 408 || size(XR, 2) ~= 200 || size(XL, 2) ~= 200 || size(XF, 2) ~= 200
            error('step0b:cols', 'Unexpected map widths seed %d', seed);
        end
        deaths = X408(:, 401:408);
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
                fit = ridge_closed(X, y);
                fprintf('  %s %s n=%d lambda=%.4g NRMSE=%.4f R2=%.4f\n', ...
                    tag, fam, nFeat, fit.lambda, fit.test_nrmse, fit.test_r2);
                rows(end+1, :) = {'mg', seed, 'driven', fam, tag, nFeat, ...
                    fit.lambda, fit.test_nrmse, fit.test_r2}; %#ok<AGROW>
            end
        end
    end

    for s = 1:numel(seeds)
        seed = seeds(s);
        folder = fullfile(benchDir, 'results', sprintf('mg_brownian_seed%d', seed));
        csvPath = fullfile(folder, 'voxels.csv');
        if ~isfile(csvPath)
            fprintf('SKIP missing %s\n', csvPath);
            skipped{end+1} = csvPath; %#ok<AGROW>
            continue
        end
        fprintf('\n=== MG brownian seed %d ===\n', seed);
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
            fit = ridge_closed(Dp, y);
            fprintf('  %s Fden n=%d lambda=%.4g NRMSE=%.4f\n', ...
                tag, nFeat, fit.lambda, fit.test_nrmse);
            rows(end+1, :) = {'mg', seed, 'brownian', 'Fden', tag, nFeat, ...
                fit.lambda, fit.test_nrmse, fit.test_r2}; %#ok<AGROW>
        end
    end

    T = cell2table(rows, 'VariableNames', {
        'task', 'seed', 'arm', 'family', 'pool', 'n_features', ...
        'lambda', 'test_nrmse', 'test_r2'});
    if ~exist(outDir, 'dir')
        mkdir(outDir);
    end
    outCsv = fullfile(outDir, 'step0b_mg_pool.csv');
    writetable(T, outCsv);
    fprintf('\nWrote %s (%d rows)\n', outCsv, height(T));
    if ~isempty(skipped)
        fprintf('Skipped %d missing path(s)\n', numel(skipped));
    end
    check_mg_sanity(T);
    report_field_vs_driven(T);
end

function check_mg_sanity(T)
    expDriven = [0.3306, 0.2152, 0.2347];
    seeds = [101, 202, 303];
    for i = 1:3
        got = pick(T, seeds(i), 'driven', 'F408', 'P0');
        if isempty(got)
            error('step0b:sanity', 'STOP missing driven F408 P0 seed %d', seeds(i));
        end
        if abs(got - expDriven(i)) > 1e-3
            error('step0b:sanity', ...
                'STOP F408 P0 seed %d NRMSE=%.6f expected %.4f', seeds(i), got, expDriven(i));
        end
        fld = pick(T, seeds(i), 'driven', 'Ffield', 'P0');
        if abs(fld - 0.0229) > 1e-3
            error('step0b:sanity', ...
                'STOP Ffield P0 seed %d NRMSE=%.6f expected 0.0229', seeds(i), fld);
        end
        den = pick(T, seeds(i), 'brownian', 'Fden', 'P0');
        if abs(den - 1.0418) > 5e-3
            error('step0b:sanity', ...
                'STOP Fden P0 seed %d NRMSE=%.6f expected ~1.0418', seeds(i), den);
        end
    end
    fprintf('SANITY PASS  F408 0.3306/0.2152/0.2347  field 0.0229  Den ~1.0418\n');
end

function report_field_vs_driven(T)
    pools = {'P0', 'P1', 'P2', 'P3'};
    allFieldWins = true;
    for i = 1:numel(pools)
        d = mean(T.test_nrmse(strcmp(T.arm, 'driven') & strcmp(T.family, 'F408') & strcmp(T.pool, pools{i})));
        f = mean(T.test_nrmse(strcmp(T.arm, 'driven') & strcmp(T.family, 'Ffield') & strcmp(T.pool, pools{i})));
        fprintf('  %s F408=%.4f  Ffield=%.4f  field_wins=%d\n', pools{i}, d, f, f < d);
        if f >= d
            allFieldWins = false;
            if ~strcmp(pools{i}, 'P3')
                error('step0b:swap', ...
                    'STOP driven F408 mean %.6f <= field %.6f at %s — check family swap', d, f, pools{i});
            end
        end
    end
    if allFieldWins
        fprintf('MG: field still wins at every pool: yes\n');
    else
        fprintf('MG: field still wins at every pool: no\n');
    end
end

function v = pick(T, seed, arm, family, pool)
    hit = T.seed == seed & strcmp(T.arm, arm) & strcmp(T.family, family) & strcmp(T.pool, pool);
    if sum(hit) == 0
        v = [];
        return
    end
    if sum(hit) ~= 1
        error('step0b:pick', 'expected 1 row for %s %s %s %d', arm, family, pool, seed);
    end
    v = T.test_nrmse(hit);
end
