function T = step2_rate_scout(varargin)
%STEP2_RATE_SCOUT  Closed ridge + occupancy for Sweep Step 2.
%   Claim seed 111 is Narma10b on disk (no BSim rerun).
%   New rates read SweepS2 results. ridge_closed unchanged. K is not a knob.

    here = fileparts(mfilename('fullpath'));
    addpath(here);
    examples = fileparts(fileparts(here));
    narmaDir = fullfile(examples, 'BSimReservoirPlanNarma10b');
    sweepDir = fullfile(examples, 'BSimReservoirPlanSweepS2');
    outDir = fullfile(sweepDir, 'results');
    targetCsv = fullfile(sweepDir, 'narma10_target.csv');
    if ~isfile(targetCsv)
        targetCsv = fullfile(narmaDir, 'narma10_target.csv');
    end

    p = inputParser;
    addParameter(p, 'NarmaDir', narmaDir);
    addParameter(p, 'SweepDir', sweepDir);
    addParameter(p, 'OutDir', outDir);
    addParameter(p, 'Target', targetCsv);
    parse(p, varargin{:});

    target = readtable(p.Results.Target, 'Delimiter', ';', 'FileType', 'text');
    y = target.y_next(:);
    u = target.u(:);
    B = baselines(y, 'Label', 'NARMA-10 SweepS2');
    fprintf('  persistence=%.4f  intercept=%.4f\n', ...
        B.persistence_nrmse, B.intercept_nrmse);

    specs = {
        'claim', 128000000, fullfile(p.Results.NarmaDir, 'results', 'narma10b_driven_seed111'), ...
                            fullfile(p.Results.NarmaDir, 'results', 'narma10b_brownian_seed111')
        'low',   64000000,  fullfile(p.Results.SweepDir, 'results', 's2_low_driven_seed111'), ...
                            fullfile(p.Results.SweepDir, 'results', 's2_low_brownian_seed111')
        'high',  256000000, fullfile(p.Results.SweepDir, 'results', 's2_high_driven_seed111'), ...
                            fullfile(p.Results.SweepDir, 'results', 's2_high_brownian_seed111')
        };

    rows = {};
    for i = 1:size(specs, 1)
        tag = specs{i, 1};
        rate = specs{i, 2};
        drivenDir = specs{i, 3};
        brownDir = specs{i, 4};
        if ~isfile(fullfile(drivenDir, 'voxels.csv'))
            error('step2:missing', 'STOP missing %s', fullfile(drivenDir, 'voxels.csv'));
        end
        if ~isfile(fullfile(brownDir, 'voxels.csv'))
            error('step2:missing', 'STOP missing %s', fullfile(brownDir, 'voxels.csv'));
        end
        fprintf('\n=== rate %s %.0f seed 111 ===\n', tag, rate);
        Td = load_voxels(drivenDir);
        X408 = window_features(Td, 'biology');
        XF = window_features(Td, 'field');
        occ = occupancy_from_table(Td, u);
        fprintf('  occupancy %s  r(mean_R,u)=%.4f  mean_R=%.4f  frac>0.5=%.4f  r(AHL,u)=%.4f\n', ...
            occ.label, occ.r_meanR_u, occ.mean_R, occ.frac_R_gt_0_5, occ.r_meanAHL_u);

        fitD = ridge_closed(X408, y);
        fitF = ridge_closed(XF, y);
        fprintf('  driven F408 n=%d lambda=%.4g NRMSE=%.4f\n', size(X408,2), fitD.lambda, fitD.test_nrmse);
        fprintf('  field  n=%d lambda=%.4g NRMSE=%.4f\n', size(XF,2), fitF.lambda, fitF.test_nrmse);

        Tb = load_voxels(brownDir);
        XD = window_features(Tb, 'brownian');
        fitB = ridge_closed(XD, y);
        fprintf('  brownian Den n=%d lambda=%.4g NRMSE=%.4f\n', size(XD,2), fitB.lambda, fitB.test_nrmse);

        rows(end+1, :) = {tag, rate, 111, 'driven', size(X408,2), fitD.lambda, ...
            fitD.test_nrmse, fitD.test_r2, occ.r_meanR_u, occ.mean_R, ...
            occ.mean_AHL, occ.frac_R_gt_0_5, occ.label}; %#ok<AGROW>
        rows(end+1, :) = {tag, rate, 111, 'field', size(XF,2), fitF.lambda, ...
            fitF.test_nrmse, fitF.test_r2, occ.r_meanR_u, occ.mean_R, ...
            occ.mean_AHL, occ.frac_R_gt_0_5, occ.label}; %#ok<AGROW>
        rows(end+1, :) = {tag, rate, 111, 'brownian', size(XD,2), fitB.lambda, ...
            fitB.test_nrmse, fitB.test_r2, NaN, NaN, NaN, NaN, 'n/a'}; %#ok<AGROW>
    end

    T = cell2table(rows, 'VariableNames', {
        'rate_tag', 'rate', 'seed', 'arm', 'n_features', 'lambda', 'test_nrmse', 'test_r2', ...
        'r_meanR_u', 'mean_R', 'mean_AHL', 'frac_R_gt_0_5', 'occupancy'});
    if ~exist(p.Results.OutDir, 'dir')
        mkdir(p.Results.OutDir);
    end
    outCsv = fullfile(p.Results.OutDir, 'step2_rate_scout.csv');
    writetable(T, outCsv);
    fprintf('\nWrote %s (%d rows)\n', outCsv, height(T));
    check_claim_sanity(T);
end

function occ = occupancy_from_table(T, u)
    rNames = T.Properties.VariableNames(startsWith(T.Properties.VariableNames, 'Receiver_R_'));
    aNames = T.Properties.VariableNames(startsWith(T.Properties.VariableNames, 'AHL_uM_'));
    windows = unique(T.Window, 'stable');
    nW = numel(windows);
    meanR = zeros(nW, 1);
    meanA = zeros(nW, 1);
    frac = zeros(nW, 1);
    for i = 1:nW
        block = T(T.Window == windows(i), :);
        R = mean(block{:, rNames}, 1, 'omitnan');
        A = mean(block{:, aNames}, 1, 'omitnan');
        meanR(i) = mean(R, 'omitnan');
        meanA(i) = mean(A, 'omitnan');
        frac(i) = mean(R > 0.5);
    end
    if numel(u) ~= nW
        error('step2:u', 'u length %d vs windows %d', numel(u), nW);
    end
    occ.r_meanR_u = pearson(meanR, u);
    occ.r_meanAHL_u = pearson(meanA, u);
    occ.mean_R = mean(meanR);
    occ.mean_AHL = mean(meanA);
    occ.frac_R_gt_0_5 = mean(frac);
    if abs(occ.r_meanR_u) < 0.5 || occ.mean_R < 0.05
        occ.label = 'dead';
    elseif occ.mean_R > 0.8 || occ.frac_R_gt_0_5 > 0.8
        occ.label = 'saturated';
    else
        occ.label = 'alive';
    end
end

function r = pearson(a, b)
    a = a(:); b = b(:);
    a = a - mean(a); b = b - mean(b);
    den = sqrt(sum(a.^2) * sum(b.^2));
    if den < 1e-15
        r = 0;
    else
        r = sum(a .* b) / den;
    end
end

function check_claim_sanity(T)
    d = pick(T, 'claim', 'driven');
    f = pick(T, 'claim', 'field');
    b = pick(T, 'claim', 'brownian');
    if abs(d - 0.9312) > 1e-3 || abs(f - 1.0289) > 1e-3 || abs(b - 1.1625) > 1e-3
        error('step2:sanity', ...
            'STOP claim seed 111 driven=%.6f field=%.6f brownian=%.6f', d, f, b);
    end
    hit = strcmp(T.rate_tag, 'claim') & strcmp(T.arm, 'driven');
    rR = T.r_meanR_u(hit);
    mR = T.mean_R(hit);
    if abs(rR - 0.8628) > 5e-3 || abs(mR - 0.1808) > 5e-3
        error('step2:sanity', ...
            'STOP claim occupancy r=%.6f mean_R=%.6f expected ~0.8628 / 0.1808', rR, mR);
    end
    fprintf('SANITY PASS  claim 111 driven=%.4f field=%.4f brownian=%.4f  r=%.4f mean_R=%.4f\n', ...
        d, f, b, rR, mR);
end

function v = pick(T, tag, arm)
    hit = strcmp(T.rate_tag, tag) & strcmp(T.arm, arm) & T.seed == 111;
    if sum(hit) ~= 1
        error('step2:pick', 'expected 1 row %s %s', tag, arm);
    end
    v = T.test_nrmse(hit);
end
