function S = kernel_rank(X, varargin)
%KERNEL_RANK  Roy–Vetterli effective rank (Wringe §8.2 / Dale adaptation).
%   S = KERNEL_RANK(X)
%
%   X is S_windows × N_channels (rows = windows, columns = channels).
%   Property measure on this matrix: all rows, no washout drop, no lambda,
%   no target y. Columns are z-scored with this matrix's mean / population
%   std (zero-variance columns stay 0). SVD of the standardised matrix.
%
%   Primary number: effective rank (Roy & Vetterli 2007)
%     p_i = σ_i / sum(σ),  R_eff = exp(-sum p_i ln p_i)
%   No singular-value cutoff for R_eff.
%
%   Also integer ranks at frozen θ ∈ {0.01, 0.05} of σ_max
%   (count σ_i > θ σ_max). Do not retune θ.
%
%   Dale/Wringe KR claim needs S >= N. If S < N, kr_claim is false and
%   ceiling = min(S,N) is printed. Official 408-D on 200 windows is a
%   ceiling, not a kernel-rank claim.
%
%   Name-value: 'Theta' (default [0.01 0.05]), 'Name' (printed label).

    p = inputParser;
    addParameter(p, 'Theta', [0.01, 0.05]);
    addParameter(p, 'Name', "");
    parse(p, varargin{:});
    theta = p.Results.Theta(:)';
    label = char(p.Results.Name);

    if ndims(X) > 2 %#ok<ISMAT>
        error('kernel_rank:shape', 'X must be a 2-D S-by-N matrix');
    end
    nWin = size(X, 1);
    nFeat = size(X, 2);
    if nWin < 2 || nFeat < 1
        error('kernel_rank:shape', 'X is %d x %d', nWin, nFeat);
    end

    mu = mean(X, 1);
    sg = std(X, 1, 1);  % population std, matches numpy default
    sg(sg < 1e-12) = 1;
    Z = (X - mu) ./ sg;

    sigma = svd(Z, 'econ');
    sigma = sigma(:);
    sigma(sigma < 0) = 0;
    total = sum(sigma);
    if total > 0
        pr = sigma / total;
        entropy = 0;
        for i = 1:numel(pr)
            if pr(i) > 0
                entropy = entropy + pr(i) * log(pr(i));
            end
        end
        rEff = exp(-entropy);
    else
        pr = zeros(size(sigma));
        rEff = 0;
    end

    sigMax = 0;
    if ~isempty(sigma)
        sigMax = sigma(1);
    end
    rankTheta = zeros(size(theta));
    for t = 1:numel(theta)
        if sigMax > 0
            rankTheta(t) = sum(sigma > theta(t) * sigMax);
        else
            rankTheta(t) = 0;
        end
    end

    S.S = nWin;
    S.N = nFeat;
    S.ceiling = min(nWin, nFeat);
    S.kr_claim = nWin >= nFeat;
    S.r_eff = rEff;
    S.theta = theta;
    S.rank_theta = rankTheta;
    S.rank_1pct = rank_at(rankTheta, theta, 0.01);
    S.rank_5pct = rank_at(rankTheta, theta, 0.05);
    S.sigma = sigma;
    S.p = pr;
    S.n_sigma = numel(sigma);
    S.n_zero_var = sum(std(X, 1, 1) < 1e-12);

    if strlength(label) == 0
        label = sprintf('%d x %d', nWin, nFeat);
    end
    claim = 'KR diagnostic (S >= N)';
    if ~S.kr_claim
        claim = sprintf('NOT a KR claim (S=%d < N=%d; ceiling %d)', ...
            nWin, nFeat, S.ceiling);
    end
    fprintf(['kernel_rank %s  S=%d N=%d  R_eff=%.3f  rank@1%%=%d  rank@5%%=%d\n' ...
             '  %s\n'], ...
        label, nWin, nFeat, rEff, S.rank_1pct, S.rank_5pct, claim);
end

function v = rank_at(rankTheta, theta, target)
    v = NaN;
    for i = 1:numel(theta)
        if abs(theta(i) - target) < 1e-15
            v = rankTheta(i);
            return
        end
    end
end
