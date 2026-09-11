function T = step1_layout_scout(varargin)
%STEP1_LAYOUT_SCOUT  Closed ridge + occupancy for Sweep Step 1.
%   Claim seed 111 is Narma10b on disk (no BSim rerun).
%   New layouts read SweepS1 results. ridge_closed unchanged.

    here = fileparts(mfilename('fullpath'));
    addpath(here);
    examples = fileparts(fileparts(here));
    narmaDir = fullfile(examples, 'BSimReservoirPlanNarma10b');
    sweepDir = fullfile(examples, 'BSimReservoirPlanSweepS1');
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
    B = baselines(y, 'Label', 'NARMA-10 SweepS1');
    fprintf('  persistence=%.4f  intercept=%.4f\n', ...
        B.persistence_nrmse, B.intercept_nrmse);

    specs = {
        'claim',     fullfile(p.Results.NarmaDir, 'results', 'narma10b_driven_seed111'), ...
                     fullfile(p.Results.NarmaDir, 'results', 'narma10b_brownian_seed111')
        'offcentre', fullfile(p.Results.SweepDir, 'results', 's1_offcentre_driven_seed111'), ...
                     fullfile(p.Results.SweepDir, 'results', 's1_offcentre_brownian_seed111')
        'wall',      fullfile(p.Results.SweepDir, 'results', 's1_wall_driven_seed111'), ...
                     fullfile(p.Results.SweepDir, 'results', 's1_wall_brownian_seed111')
        };

    rows = {};
    for i = 1:size(specs, 1)
        layout = specs{i, 1};
        drivenDir = specs{i, 2};
        brownDir = specs{i, 3};
        if ~isfile(fullfile(drivenDir, 'voxels.csv'))
            error('step1:missing', 'STOP missing %s', fullfile(drivenDir, 'voxels.csv'));
        end
        if ~isfile(fullfile(brownDir, 'voxels.csv'))
            error('step1:missing', 'STOP missing %s', fullfile(brownDir, 'voxels.csv'));
        end
        fprintf('\n=== layout %s seed 111 ===\n', layout);
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

        rows(end+1, :) = {layout, 111, 'driven', size(X408,2), fitD.lambda, ...
            fitD.test_nrmse, fitD.test_r2, occ.r_meanR_u, occ.mean_R, ...
            occ.mean_AHL, occ.frac_R_gt_0_5, occ.label}; %#ok<AGROW>
        rows(end+1, :) = {layout, 111, 'field', size(XF,2), fitF.lambda, ...
            fitF.test_nrmse, fitF.test_r2, occ.r_meanR_u, occ.mean_R, ...
            occ.mean_AHL, occ.frac_R_gt_0_5, occ.label}; %#ok<AGROW>
        rows(end+1, :) = {layout, 111, 'brownian', size(XD,2), fitB.lambda, ...
            fitB.test_nrmse, fitB.test_r2, NaN, NaN, NaN, NaN, 'n/a'}; %#ok<AGROW>
    end

    T = cell2table(rows, 'VariableNames', {
        'layout', 'seed', 'arm', 'n_features', 'lambda', 'test_nrmse', 'test_r2', ...
        'r_meanR_u', 'mean_R', 'mean_AHL', 'frac_R_gt_0_5', 'occupancy'});
    if ~exist(p.Results.OutDir, 'dir')
        mkdir(p.Results.OutDir);
    end
    outCsv = fullfile(p.Results.OutDir, 'step1_layout_scout.csv');
    writetable(T, outCsv);
    fprintf('\nWrote %s (%d rows)\n', outCsv, height(T));
    check_claim_sanity(T);
end

function occ = occupancy_from_table(T, u)
    rNames = T.Properties.VariableNames(startsWith(T.Properties.VariableNames, 'Receiver_R_'));
    lNames = T.Properties.VariableNames(startsWith(T.Properties.VariableNames, 'Lum_Mean_'));
    aNames = T.Properties.VariableNames(startsWith(T.Properties.VariableNames, 'AHL_uM_'));
    windows = unique(T.Window, 'stable');
    nW = numel(windows);
    meanR = zeros(nW, 1);
    meanL = zeros(nW, 1);
    meanA = zeros(nW, 1);
    frac = zeros(nW, 1);
    for i = 1:nW
        block = T(T.Window == windows(i), :);
        R = mean(block{:, rNames}, 1, 'omitnan');
        L = mean(block{:, lNames}, 1, 'omitnan');
        A = mean(block{:, aNames}, 1, 'omitnan');
        meanR(i) = mean(R, 'omitnan');
        meanL(i) = mean(L, 'omitnan');
        meanA(i) = mean(A, 'omitnan');
        frac(i) = mean(R > 0.5);
    end
    if numel(u) ~= nW
        error('step1:u', 'u length %d vs windows %d', numel(u), nW);
    end
    occ.r_meanR_u = pearson(meanR, u);
    occ.r_meanAHL_u = pearson(meanA, u);
    occ.mean_R = mean(meanR);
    occ.mean_L = mean(meanL);
    occ.mean_AHL = mean(meanA);
    occ.frac_R_gt_0_5 = mean(frac);
    if abs(occ.r_meanR_u) < 0.5 || occ.mean_R < 0.05
        occ.label = 'dead';
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
        error('step1:sanity', ...
            'STOP claim seed 111 driven=%.6f field=%.6f brownian=%.6f', d, f, b);
    end
    fprintf('SANITY PASS  claim 111 driven=%.4f field=%.4f brownian=%.4f\n', d, f, b);
end

function v = pick(T, layout, arm)
    hit = strcmp(T.layout, layout) & strcmp(T.arm, arm) & T.seed == 111;
    if sum(hit) ~= 1
        error('step1:pick', 'expected 1 row %s %s', layout, arm);
    end
    v = T.test_nrmse(hit);
end
