function S = memory_capacity(X, u, varargin)
%MEMORY_CAPACITY  Linear random-signal recall MC_k for k=1..20.
%   S = MEMORY_CAPACITY(X, u)
%   Same closed-slice ridge as ridge_closed. Reconstructs u[n-k].
%   MC = sum max(0, test R^2_k) for k=1..20. k=0 is printed, not added.
%
%   Independent lambda per delay. Do not put voxel AHL in X for biology.

    p = inputParser;
    addParameter(p, 'Kmax', 20);
    parse(p, varargin{:});
    kmax = p.Results.Kmax;

    u = u(:);
    if size(X, 1) ~= numel(u)
        error('memory_capacity:length', 'X has %d rows, u has %d', size(X, 1), numel(u));
    end
    n = numel(u);
    idx = (0:n-1)';

    delays = struct('k', {}, 'lambda', {}, 'test_nrmse', {}, 'test_r2', {}, 'contribution', {});
    mc = 0;
    for k = 1:kmax
        target = zeros(n, 1);
        ok = idx >= k;
        target(ok) = u(idx(ok) - k + 1);
        fit = ridge_closed(X, target);
        r2 = fit.test_r2;
        contrib = 0;
        if isfinite(r2)
            contrib = max(0, r2);
        end
        mc = mc + contrib;
        delays(k).k = k;
        delays(k).lambda = fit.lambda;
        delays(k).test_nrmse = fit.test_nrmse;
        delays(k).test_r2 = r2;
        delays(k).contribution = contrib;
        fprintf('  k=%2d  R^2=%.4f  NRMSE=%.4f  lambda=%.4g\n', ...
            k, r2, fit.test_nrmse, fit.lambda);
    end

    k0 = ridge_closed(X, u);
    S.mc = mc;
    S.kmax = kmax;
    S.delays = delays;
    S.k0_test_r2 = k0.test_r2;
    S.k0_test_nrmse = k0.test_nrmse;
    S.k0_lambda = k0.lambda;
    fprintf('Linear MC (k=1..%d) = %.4f   k=0 R^2 = %.4f (not added)\n', ...
        kmax, S.mc, S.k0_test_r2);
end
