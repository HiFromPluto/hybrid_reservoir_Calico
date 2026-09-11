function S = horizon_narma(X, y_next, h, varargin)
%HORIZON_NARMA  Teacher-forced h-step NARMA readout. Not closed-loop free-run.
%   S = HORIZON_NARMA(X, y_next, h)
%
%   Window n predicts y[n+h] = y_next[n+h-1]. u still drove the dish.
%   Rows with n+h-1 > 199 are dropped. Train windows 40..149 stay valid
%   for h in {1,5,10}. Test is 150..(199-h+1); n_test = 50 / 46 / 41.
%   Independent lambda. Same closed val-slice as ridge_closed.

    if nargin < 3 || isempty(h)
        h = 1;
    end
    y_next = y_next(:);
    n = size(X, 1);
    if n ~= numel(y_next)
        error('horizon_narma:length', 'X has %d rows, y_next has %d', n, numel(y_next));
    end
    if h < 1 || h > n
        error('horizon_narma:h', 'h must be in 1..%d', n);
    end

    nValid = n - h + 1;
    yH = y_next(h:h+nValid-1);
    XH = X(1:nValid, :);
    fit = ridge_closed(XH, yH);

    S.h = h;
    S.n_valid = nValid;
    S.n_test = nValid - 40 - 110;
    S.lambda = fit.lambda;
    S.train_nrmse = fit.train_nrmse;
    S.test_nrmse = fit.test_nrmse;
    S.test_r2 = fit.test_r2;
    S.y_test = fit.y_test;
    S.y_pred = fit.y_pred;
    fprintf('Teacher-forced horizon h=%d  n_test=%d  NRMSE=%.4f  R^2=%.4f  lambda=%.4g\n', ...
        h, S.n_test, S.test_nrmse, S.test_r2, S.lambda);
end
