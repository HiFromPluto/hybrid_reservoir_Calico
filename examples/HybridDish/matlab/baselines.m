function S = baselines(y, varargin)
%BASELINES  Intercept and teacher-forced persistence NRMSE.
%   S = BASELINES(y)                 y is 200 x 1 target (y_next or x_next)
%   S = BASELINES(target_csv)
%
%   Intercept: constant = mean of train target (windows 40..149).
%   Persistence: window n predicts the next sample by the current one
%   (y[n] = y_next[n-1] / x[n] = x_next[n-1]).
%   Both evaluated on test windows 150..199. NRMSE uses test pop-std.
%
%   These baselines do not use voxels. They are Wringe comparators for the
%   frozen target, not a rewrite of any GATE Overall line.
%
%   Name-value: 'Label' (default 'NARMA-10') is printed in the fprintf.

    p = inputParser;
    addParameter(p, 'Washout', 40);
    addParameter(p, 'Train', 110);
    addParameter(p, 'Test', 50);
    addParameter(p, 'Label', 'NARMA-10');
    parse(p, varargin{:});
    washout = p.Results.Washout;
    nTrain = p.Results.Train;
    nTest = p.Results.Test;
    label = char(p.Results.Label);

    if ischar(y) || isstring(y)
        T = readtable(char(y), 'Delimiter', ';', 'FileType', 'text');
        if ismember('x_next', T.Properties.VariableNames)
            y = T.x_next;
        else
            y = T.y_next;
        end
    end
    y = y(:);
    n = numel(y);
    if n < washout + nTrain + nTest
        error('baselines:length', 'y_next is too short (%d)', n);
    end

    train = y(washout + (1:nTrain));
    test = y(washout + nTrain + (1:nTest));
    trainMean = mean(train);

    interceptPred = trainMean * ones(size(test));
    persistPred = y(washout + nTrain + (1:nTest) - 1);  % y_next[n-1] for window n

    S.n_test = nTest;
    S.train_mean = trainMean;
    S.test_mean = mean(test);
    S.test_std = std(test, 1);
    S.intercept_nrmse = nrmse_pop(test, interceptPred);
    S.persistence_nrmse = nrmse_pop(test, persistPred);
    S.intercept_pred = interceptPred;
    S.persistence_pred = persistPred;
    S.y_test = test;

    fprintf('%s baselines (test %d windows)\n', label, nTest);
    fprintf('  intercept   NRMSE = %.4f  (train mean %.4f, test mean %.4f)\n', ...
        S.intercept_nrmse, S.train_mean, S.test_mean);
    fprintf('  persistence NRMSE = %.4f  (teacher-forced current -> next)\n', ...
        S.persistence_nrmse);
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
