function S = ridge_closed(X, y, varargin)
%RIDGE_CLOSED  Ridge with washout / train / test and a closed validation slice.
%   Lambda is chosen on train windows 88..109 only (after washout), never on test.
%
%   Name-value:
%     'Washout'  (default 40)
%     'Train'    (default 110)
%     'Test'     (default 50)
%     'InnerVal' (default 22)
%     'Grid'     (default [1e-6 1e-4 1e-2 1 1e2 1e4 1e6])

    p = inputParser;
    addParameter(p, 'Washout', 40);
    addParameter(p, 'Train', 110);
    addParameter(p, 'Test', 50);
    addParameter(p, 'InnerVal', 22);
    addParameter(p, 'Grid', [1e-6, 1e-4, 1e-2, 1, 1e2, 1e4, 1e6]);
    parse(p, varargin{:});
    washout = p.Results.Washout;
    nTrain = p.Results.Train;
    nTest = p.Results.Test;
    innerVal = p.Results.InnerVal;
    grid = p.Results.Grid(:)';

    y = y(:);
    if size(X, 1) ~= numel(y)
        error('ridge_closed:length', 'X has %d rows, y has %d', size(X, 1), numel(y));
    end
    inferredTest = size(X, 1) - washout - nTrain;
    if inferredTest < 1
        error('ridge_closed:length', ...
            'Need washout+train+%d rows, got %d', nTest, size(X, 1));
    end
    % Horizon tasks drop the last h-1 windows; infer test length when X is shorter.
    if inferredTest ~= nTest
        nTest = inferredTest;
    end

    Xw = X(washout+1:end, :);
    yw = y(washout+1:end);

    innerEnd = nTrain - innerVal;
    Xinner = Xw(1:innerEnd, :);
    yinner = yw(1:innerEnd);
    Xval = Xw(innerEnd+1:nTrain, :);
    yval = yw(innerEnd+1:nTrain);

    [XinnerZ, mu, sg] = standardize_fit(Xinner);
    XvalZ = standardize_apply(Xval, mu, sg);

    bestScore = inf;
    bestLam = grid(end);
    valScores = zeros(size(grid));
    for i = 1:numel(grid)
        w = ridge_fit(XinnerZ, yinner, grid(i));
        pred = ridge_predict(XvalZ, w);
        valScores(i) = nrmse_pop(yval, pred);
        % Lower NRMSE wins; ties take the larger lambda.
        if valScores(i) < bestScore - 1e-15 || ...
                (abs(valScores(i) - bestScore) <= 1e-15 && grid(i) > bestLam)
            bestScore = valScores(i);
            bestLam = grid(i);
        end
    end

    Xtrain = Xw(1:nTrain, :);
    ytrain = yw(1:nTrain);
    Xtest = Xw(nTrain+1:nTrain+nTest, :);
    ytest = yw(nTrain+1:nTrain+nTest);
    [XtrainZ, mu, sg] = standardize_fit(Xtrain);
    XtestZ = standardize_apply(Xtest, mu, sg);
    weights = ridge_fit(XtrainZ, ytrain, bestLam);
    trainPred = ridge_predict(XtrainZ, weights);
    testPred = ridge_predict(XtestZ, weights);

    S.lambda = bestLam;
    S.val_scores = [grid(:), valScores(:)];
    S.train_nrmse = nrmse_pop(ytrain, trainPred);
    S.test_nrmse = nrmse_pop(ytest, testPred);
    S.test_r2 = r_squared(ytest, testPred);
    S.y_test = ytest;
    S.y_pred = testPred;
    S.weights = weights;
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

function v = nrmse_pop(ytrue, ypred)
    ytrue = ytrue(:);
    ypred = ypred(:);
    rmse = sqrt(mean((ytrue - ypred).^2));
    denom = std(ytrue, 1);
    if denom > 0
        v = rmse / denom;
    else
        v = inf;
    end
end

function v = r_squared(ytrue, ypred)
    ytrue = ytrue(:);
    ypred = ypred(:);
    ssRes = sum((ytrue - ypred).^2);
    ssTot = sum((ytrue - mean(ytrue)).^2);
    if ssTot > 0
        v = 1 - ssRes / ssTot;
    else
        v = NaN;
    end
end
