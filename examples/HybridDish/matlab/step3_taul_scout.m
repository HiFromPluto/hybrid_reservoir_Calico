function T = step3_taul_scout(varargin)
%STEP3_TAUL_SCOUT  Closed ridge + L lag diagnostics for Sweep Step 3.
%   Claim seed 111 is Narma10b on disk. Brownian reused from claim.
%   ridge_closed unchanged. tau_R / K / rate / position are not knobs.

    here = fileparts(mfilename('fullpath'));
    addpath(here);
    examples = fileparts(fileparts(here));
    narmaDir = fullfile(examples, 'BSimReservoirPlanNarma10b');
    sweepDir = fullfile(examples, 'BSimReservoirPlanSweepS3');
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
    B = baselines(y, 'Label', 'NARMA-10 SweepS3');
    fprintf('  persistence=%.4f  intercept=%.4f\n', ...
        B.persistence_nrmse, B.intercept_nrmse);

    claimBrown = fullfile(p.Results.NarmaDir, 'results', 'narma10b_brownian_seed111');
    if ~isfile(fullfile(claimBrown, 'voxels.csv'))
        error('step3:missing', 'STOP missing %s', fullfile(claimBrown, 'voxels.csv'));
    end
    Tb = load_voxels(claimBrown);
    XD = window_features(Tb, 'brownian');
    fitB = ridge_closed(XD, y);
    fprintf('  reused claim brownian Den n=%d lambda=%.4g NRMSE=%.4f\n', ...
        size(XD, 2), fitB.lambda, fitB.test_nrmse);

    specs = {
        'claim', 1500, 1/1500, 1/1500, ...
            fullfile(p.Results.NarmaDir, 'results', 'narma10b_driven_seed111')
        'fast',  500,  1/500,  1/500, ...
            fullfile(p.Results.SweepDir, 'results', 's3_fast_driven_seed111')
        'slow',  4500, 1/4500, 1/4500, ...
            fullfile(p.Results.SweepDir, 'results', 's3_slow_driven_seed111')
        };

    rows = {};
    for i = 1:size(specs, 1)
        tag = specs{i, 1};
        tauL = specs{i, 2};
        alpha = specs{i, 3};
        delta = specs{i, 4};
        drivenDir = specs{i, 5};
        if ~isfile(fullfile(drivenDir, 'voxels.csv'))
            error('step3:missing', 'STOP missing %s', fullfile(drivenDir, 'voxels.csv'));
        end
        fprintf('\n=== tau_L %s %.0f s seed 111 ===\n', tag, tauL);
        Td = load_voxels(drivenDir);
        X408 = window_features(Td, 'biology');
        XF = window_features(Td, 'field');
        XL = window_features(Td, 'lum');
        XR = window_features(Td, 'receiver');
        diag = lag_diagnostics(Td, u);
        fprintf('  R r(u)=%.4f mean_R=%.4f\n', diag.r_meanR_u, diag.mean_R);
        fprintf('  L r(u)=%.4f r(L,R)=%.4f mean_L=%.4f  L[0]=%.4f L[40]=%.4f L[199]=%.4f\n', ...
            diag.r_meanL_u, diag.r_meanL_meanR, diag.mean_L, ...
            diag.mean_L_w0, diag.mean_L_w40, diag.mean_L_w199);
        fprintf('  AHL r(u)=%.4f\n', diag.r_meanAHL_u);
        if abs(diag.r_meanR_u) < 0.5
            error('step3:receiver', ...
                'STOP RECEIVER_BROKEN %s r(mean_R,u)=%.4f — wrong constant', tag, diag.r_meanR_u);
        end

        fitD = ridge_closed(X408, y);
        fitF = ridge_closed(XF, y);
        fitL = ridge_closed(XL, y);
        fitR = ridge_closed(XR, y);
        fprintf('  driven F408 n=%d lambda=%.4g NRMSE=%.4f\n', size(X408,2), fitD.lambda, fitD.test_nrmse);
        fprintf('  field  n=%d lambda=%.4g NRMSE=%.4f\n', size(XF,2), fitF.lambda, fitF.test_nrmse);
        fprintf('  FL     n=%d lambda=%.4g NRMSE=%.4f\n', size(XL,2), fitL.lambda, fitL.test_nrmse);
        fprintf('  FR     n=%d lambda=%.4g NRMSE=%.4f\n', size(XR,2), fitR.lambda, fitR.test_nrmse);

        rows(end+1, :) = row_from(tag, tauL, alpha, delta, 'driven', fitD, size(X408,2), diag); %#ok<AGROW>
        rows(end+1, :) = row_from(tag, tauL, alpha, delta, 'field', fitF, size(XF,2), diag); %#ok<AGROW>
        rows(end+1, :) = {tag, tauL, alpha, delta, 111, 'brownian', size(XD,2), fitB.lambda, ...
            fitB.test_nrmse, fitB.test_r2, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN}; %#ok<AGROW>
        rows(end+1, :) = row_from(tag, tauL, alpha, delta, 'FL', fitL, size(XL,2), diag); %#ok<AGROW>
        rows(end+1, :) = row_from(tag, tauL, alpha, delta, 'FR', fitR, size(XR,2), diag); %#ok<AGROW>
    end

    T = cell2table(rows, 'VariableNames', {
        'tau_tag', 'tau_L', 'alpha', 'delta', 'seed', 'arm', 'n_features', ...
        'lambda', 'test_nrmse', 'test_r2', 'r_meanR_u', 'r_meanL_u', 'r_meanL_meanR', ...
        'mean_R', 'mean_L', 'mean_L_w0', 'mean_L_w40', 'mean_L_w199'});
    if ~exist(p.Results.OutDir, 'dir')
        mkdir(p.Results.OutDir);
    end
    outCsv = fullfile(p.Results.OutDir, 'step3_taul_scout.csv');
    writetable(T, outCsv);
    fprintf('\nWrote %s (%d rows)\n', outCsv, height(T));
    check_claim_sanity(T);
end

function r = row_from(tag, tauL, alpha, delta, arm, fit, nFeat, diag)
    r = {tag, tauL, alpha, delta, 111, arm, nFeat, fit.lambda, fit.test_nrmse, fit.test_r2, ...
        diag.r_meanR_u, diag.r_meanL_u, diag.r_meanL_meanR, diag.mean_R, diag.mean_L, ...
        diag.mean_L_w0, diag.mean_L_w40, diag.mean_L_w199};
end

function d = lag_diagnostics(T, u)
    rNames = T.Properties.VariableNames(startsWith(T.Properties.VariableNames, 'Receiver_R_'));
    lNames = T.Properties.VariableNames(startsWith(T.Properties.VariableNames, 'Lum_Mean_'));
    aNames = T.Properties.VariableNames(startsWith(T.Properties.VariableNames, 'AHL_uM_'));
    windows = unique(T.Window, 'stable');
    nW = numel(windows);
    meanR = zeros(nW, 1);
    meanL = zeros(nW, 1);
    meanA = zeros(nW, 1);
    for i = 1:nW
        block = T(T.Window == windows(i), :);
        R = mean(block{:, rNames}, 1, 'omitnan');
        L = mean(block{:, lNames}, 1, 'omitnan');
        A = mean(block{:, aNames}, 1, 'omitnan');
        meanR(i) = mean(R, 'omitnan');
        meanL(i) = mean(L, 'omitnan');
        meanA(i) = mean(A, 'omitnan');
    end
    if numel(u) ~= nW
        error('step3:u', 'u length %d vs windows %d', numel(u), nW);
    end
    d.r_meanR_u = pearson(meanR, u);
    d.r_meanL_u = pearson(meanL, u);
    d.r_meanL_meanR = pearson(meanL, meanR);
    d.r_meanAHL_u = pearson(meanA, u);
    d.mean_R = mean(meanR);
    d.mean_L = mean(meanL);
    d.mean_L_w0 = meanL(windows == 0);
    d.mean_L_w40 = meanL(windows == 40);
    d.mean_L_w199 = meanL(windows == 199);
    if isempty(d.mean_L_w0) || isempty(d.mean_L_w40) || isempty(d.mean_L_w199)
        error('step3:windows', 'Need mean_L at windows 0, 40, 199');
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
        error('step3:sanity', ...
            'STOP claim seed 111 driven=%.6f field=%.6f brownian=%.6f', d, f, b);
    end
    fprintf('SANITY PASS  claim 111 driven=%.4f field=%.4f brownian=%.4f\n', d, f, b);
end

function v = pick(T, tag, arm)
    hit = strcmp(T.tau_tag, tag) & strcmp(T.arm, arm) & T.seed == 111;
    if sum(hit) ~= 1
        error('step3:pick', 'expected 1 row %s %s', tag, arm);
    end
    v = T.test_nrmse(hit);
end
