function S = classify_waveform(voxel_dir, labels_csv, varargin)
%CLASSIFY_WAVEFORM  One-vs-rest ridge on sine/square/triangle labels.
%   Matches examples/BSimReservoirPlanWaveform/check_waveform.py.
%
%   One lambda for all three heads, chosen by highest macro OVR
%   validation AUC on post-washout rows 88..109. Ties take the larger
%   lambda. Do not pick lambda by NRMSE on 0/1 labels.
%
%   Chance accuracy is 1/3. Chance macro OVR AUC is 0.5.
%   Official biology features are R + L + deaths, not voxel AHL.
%
%   Name-value:
%     'Family'  'biology' (default) | 'brownian' | 'field'

    if nargin < 2 || isempty(labels_csv)
        here = fileparts(mfilename('fullpath'));
        labels_csv = fullfile(fileparts(fileparts(here)), ...
            'BSimReservoirPlanWaveform', 'waveform_labels.csv');
        if ~isfile(labels_csv)
            labels_csv = fullfile(fileparts(here), 'input', 'waveform_labels.csv');
        end
    end
    p = inputParser;
    addParameter(p, 'Family', 'biology');
    parse(p, varargin{:});

    T = load_voxels(voxel_dir);
    [X, names] = window_features(T, p.Results.Family);
    L = readtable(labels_csv, 'Delimiter', ';', 'FileType', 'text');
    y = L.y(:);
    if numel(y) ~= size(X, 1)
        error('classify_waveform:windows', 'Label/window length mismatch');
    end

    washout = 40;
    nTrain = 110;
    nTest = 50;
    innerVal = 22;
    grid = [1e-6, 1e-4, 1e-2, 1, 1e2, 1e4, 1e6];
    nClass = 3;
    classNames = {'sine', 'square', 'triangle'};

    Xw = X(washout+1:end, :);
    yw = y(washout+1:end);
    if size(Xw, 1) ~= nTrain + nTest
        error('classify_waveform:split', 'expected 160 post-washout rows');
    end

    innerEnd = nTrain - innerVal;
    Xinner = Xw(1:innerEnd, :);
    yinner = yw(1:innerEnd);
    Xval = Xw(innerEnd+1:nTrain, :);
    yval = yw(innerEnd+1:nTrain);
    [XinnerZ, muInner, sgInner] = standardize_fit(Xinner);
    XvalZ = standardize_apply(Xval, muInner, sgInner);

    bestAuc = -inf;
    bestLam = grid(end);
    valGrid = zeros(numel(grid), 1);
    for i = 1:numel(grid)
        scoresVal = ovr_scores(XinnerZ, yinner, XvalZ, grid(i), nClass);
        [macro, ~] = macro_ovr_auc(yval, scoresVal, nClass);
        valGrid(i) = macro;
        if macro > bestAuc + 1e-15 || ...
                (abs(macro - bestAuc) <= 1e-15 && grid(i) > bestLam)
            bestAuc = macro;
            bestLam = grid(i);
        end
    end

    Xtrain = Xw(1:nTrain, :);
    ytrain = yw(1:nTrain);
    Xtest = Xw(nTrain+1:nTrain+nTest, :);
    ytest = yw(nTrain+1:nTrain+nTest);
    [XtrainZ, mu, sg] = standardize_fit(Xtrain);
    XtestZ = standardize_apply(Xtest, mu, sg);
    scores = ovr_scores(XtrainZ, ytrain, XtestZ, bestLam, nClass);
    [macroAuc, aucPer] = macro_ovr_auc(ytest, scores, nClass);
    [~, pred] = max(scores, [], 2);
    pred = pred - 1;  % classes 0/1/2
    acc = mean(pred == ytest);
    C = confusion_counts(ytest, pred, nClass);
    prf = prf_from_confusion(C);

    S.family = p.Results.Family;
    S.n_features = numel(names);
    S.lambda = bestLam;
    S.val_macro_ovr_auc = valGrid;
    S.val_grid = grid(:);
    S.accuracy = acc;
    S.macro_auc = macroAuc;
    S.auc_per_class = aucPer(:);
    S.classes = (0:nClass-1)';
    S.class_names = classNames;
    S.y_test = ytest;
    S.y_pred = pred;
    S.scores = scores;
    S.confusion = C;
    S.precision = prf.precision;
    S.recall = prf.recall;
    S.f1 = prf.f1;
    S.support = prf.support;
    S.macro_f1 = prf.macro_f1;
    S.chance_accuracy = 1/3;
    S.chance_auc = 0.5;

    fprintf('HybridDish waveform (%s, %d features, lambda=%.4g)\n', ...
        S.family, S.n_features, S.lambda);
    fprintf('  accuracy   = %.4f  (chance 0.333)\n', S.accuracy);
    fprintf('  macro AUC  = %.4f  (chance 0.500)\n', S.macro_auc);
    fprintf('  macro F1   = %.4f\n', S.macro_f1);
    fprintf('  confusion  = [%s; %s; %s]\n', ...
        mat2str(C(1,:)), mat2str(C(2,:)), mat2str(C(3,:)));
    for k = 1:nClass
        fprintf('  %-8s P=%.3f R=%.3f F1=%.3f support=%d AUC=%.3f\n', ...
            classNames{k}, S.precision(k), S.recall(k), S.f1(k), ...
            S.support(k), S.auc_per_class(k));
    end
end

function scores = ovr_scores(X, y, Xpred, lam, nClass)
    scores = zeros(size(Xpred, 1), nClass);
    for k = 0:nClass-1
        w = ridge_fit(X, double(y == k), lam);
        scores(:, k+1) = ridge_predict(Xpred, w);
    end
end

function [macro, per] = macro_ovr_auc(y, scores, nClass)
    per = zeros(nClass, 1);
    for k = 0:nClass-1
        per(k+1) = mann_whitney_auc(y == k, scores(:, k+1));
    end
    macro = mean(per);
end

function a = mann_whitney_auc(ytrue, scores)
% Mann-Whitney / Wilcoxon on scores vs {0,1}; average ranks for ties.
    ytrue = ytrue(:) ~= 0;
    scores = scores(:);
    nPos = sum(ytrue);
    nNeg = sum(~ytrue);
    if nPos == 0 || nNeg == 0
        a = NaN;
        return
    end
    [sorted, order] = sort(scores, 'ascend');
    ranks = zeros(size(scores));
    i = 1;
    n = numel(scores);
    while i <= n
        j = i;
        while j < n && sorted(j+1) == sorted(i)
            j = j + 1;
        end
        ranks(order(i:j)) = 0.5 * (i + j);
        i = j + 1;
    end
    posRankSum = sum(ranks(ytrue));
    a = (posRankSum - nPos * (nPos + 1) / 2) / (nPos * nNeg);
end

function C = confusion_counts(ytrue, ypred, nClass)
    C = zeros(nClass, nClass);
    for i = 1:numel(ytrue)
        C(ytrue(i)+1, ypred(i)+1) = C(ytrue(i)+1, ypred(i)+1) + 1;
    end
end

function prf = prf_from_confusion(C)
    nClass = size(C, 1);
    precision = zeros(nClass, 1);
    recall = zeros(nClass, 1);
    f1 = zeros(nClass, 1);
    support = sum(C, 2);
    for k = 1:nClass
        tp = C(k, k);
        fp = sum(C(:, k)) - tp;
        fn = sum(C(k, :)) - tp;
        if tp + fp > 0
            precision(k) = tp / (tp + fp);
        else
            precision(k) = 0;
        end
        if tp + fn > 0
            recall(k) = tp / (tp + fn);
        else
            recall(k) = 0;
        end
        if precision(k) + recall(k) > 0
            f1(k) = 2 * precision(k) * recall(k) / (precision(k) + recall(k));
        else
            f1(k) = 0;
        end
    end
    prf.precision = precision;
    prf.recall = recall;
    prf.f1 = f1;
    prf.support = support;
    prf.macro_f1 = mean(f1);
end

function w = ridge_fit(X, y, lam)
    n = size(X, 2);
    Xb = [ones(size(X, 1), 1), X];
    G = Xb' * Xb;
    G(2:end, 2:end) = G(2:end, 2:end) + lam * eye(n);
    rhs = Xb' * y(:);
    w = G \ rhs;
end

function yhat = ridge_predict(X, w)
    Xb = [ones(size(X, 1), 1), X];
    yhat = Xb * w;
end

function [Z, mu, sg] = standardize_fit(X)
    mu = mean(X, 1);
    sg = std(X, 1, 1);  % population std, matches numpy default
    sg(sg < 1e-12) = 1;
    Z = (X - mu) ./ sg;
end

function Z = standardize_apply(X, mu, sg)
    Z = (X - mu) ./ sg;
end
