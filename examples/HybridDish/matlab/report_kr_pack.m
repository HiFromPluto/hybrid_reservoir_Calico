function R = report_kr_pack(varargin)
%REPORT_KR_PACK  Wringe pack 4 kernel-rank diagnostic on Narma10b voxels.
%
%   From the HybridDish folder:
%     addpath('matlab')
%     R = report_kr_pack
%
%   One driven seed (MATLAB reproduction target):
%     R = report_kr_pack('Only', "driven_seed222")
%
%   Frozen Narma10b voxels. Does not run BSim. No lambda, no target y.
%   All 200 windows. Primary families are 20x10 maps (S = N = 200).
%   Official 408-D is printed as a ceiling, not a KR claim.
%   Python report_kr_pack.py writes WRINGE_KR_PACK.md.

    here = fileparts(mfilename('fullpath'));
    dish = fileparts(here);
    examples = fileparts(dish);
    narmaDir = fullfile(examples, 'BSimReservoirPlanNarma10b');

    p = inputParser;
    addParameter(p, 'NarmaDir', narmaDir);
    addParameter(p, 'Only', "");
    parse(p, varargin{:});

    addpath(here);
    arms = {'driven', 'brownian', 'silent'};
    seeds = [111, 222, 333];
    R.kr = struct();

    for a = 1:numel(arms)
        arm = arms{a};
        for s = 1:numel(seeds)
            tag = sprintf('%s_seed%d', arm, seeds(s));
            if strlength(p.Results.Only) > 0 && ~strcmp(tag, char(p.Results.Only))
                continue
            end
            folder = fullfile(p.Results.NarmaDir, 'results', ['narma10b_' tag]);
            fprintf('\n=== KR Narma10b %s ===\n', tag);
            R.kr.(matlab.lang.makeValidName(tag)) = analyse_one(folder, arm);
        end
    end
end

function S = analyse_one(folder, arm)
    T = load_voxels(folder);
    S.folder = folder;
    S.arm = arm;
    if strcmp(arm, 'brownian')
        [X, names] = window_features(T, 'brownian');
        S.den = kernel_rank(X, 'Name', sprintf('Den_* n=%d', numel(names)));
        return
    end
    [Xl, nl] = window_features(T, 'lum');
    S.lum = kernel_rank(Xl, 'Name', sprintf('Lum_Mean_* n=%d', numel(nl)));
    [Xr, nr] = window_features(T, 'receiver');
    S.receiver = kernel_rank(Xr, 'Name', sprintf('Receiver_R_* n=%d', numel(nr)));
    [Xb, nb] = window_features(T, 'biology');
    S.biology408 = kernel_rank(Xb, 'Name', sprintf('official 408 n=%d', numel(nb)));
    if strcmp(arm, 'driven')
        [Xf, nf] = window_features(T, 'field');
        S.field = kernel_rank(Xf, 'Name', sprintf('AHL_uM_* n=%d', numel(nf)));
    end
end
