function S = classify_waveform_features(X, y, varargin)
%CLASSIFY_WAVEFORM_FEATURES  OVR ridge on already-built window features.
%   Same rules as classify_waveform.m. Do not pick lambda by NRMSE.
%   One lambda for all three heads: highest macro OVR val AUC on
%   post-washout rows 88..109. Ties take the larger lambda.

    p = inputParser;
    addParameter(p, 'Quiet', false);
    parse(p, varargin{:});

    y = y(:);
    if numel(y) ~= size(X, 1)
        error('classify_waveform_features:windows', 'Label/window length mismatch');
    end

    washout = 40;
    nTrain = 110;
    nTest = 50;
    innerVal = 22;
    grid = [1e-6, 1e-4, 1e-2, 1, 1e2, 1e4, 1e6];
    nClass = 3;

    Xw = X(washout+1:end, :);
    yw = y(washout+1:end);
    if size(Xw, 1) ~= nTrain + nTest
        error('classify_waveform_features:split', 'expected 160 post-washout rows');
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
    [macroAuc, ~] = macro_ovr_auc(ytest, scores, nClass);
    [~, pred] = max(scores, [], 2);
    pred = pred - 1;
    acc = mean(pred == ytest);

    S.n_features = size(X, 2);
    S.lambda = bestLam;
    S.val_macro_ovr_auc = valGrid;
    S.accuracy = acc;
    S.macro_auc = macroAuc;
    S.chance_accuracy = 1/3;
    S.chance_auc = 0.5;

    if ~p.Results.Quiet
        fprintf('  n=%d lambda=%.4g AUC=%.4f acc=%.4f\n', ...
            S.n_features, S.lambda, S.macro_auc, S.accuracy);
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
    sg = std(X, 1, 1);
    sg(sg < 1e-12) = 1;
    Z = (X - mu) ./ sg;
end

function Z = standardize_apply(X, mu, sg)
    Z = (X - mu) ./ sg;
end
