function R = report_waveform_pack(varargin)
%REPORT_WAVEFORM_PACK  Wringe pack 2 driver for the nine waveform dirs.
%
%   From the HybridDish folder:
%     addpath('matlab')
%     R = report_waveform_pack
%
%   One driven seed (MATLAB reproduction target):
%     R = report_waveform_pack('Only', "driven_seed101")
%
%   Frozen waveform voxels. Does not run BSim. Val-AUC lambda, not NRMSE.
%   Python report_waveform_pack.py writes WRINGE_WAVEFORM_PACK.md.

    here = fileparts(mfilename('fullpath'));
    dish = fileparts(here);
    examples = fileparts(dish);
    waveDir = fullfile(examples, 'BSimReservoirPlanWaveform');
    labelsCsv = fullfile(waveDir, 'waveform_labels.csv');

    p = inputParser;
    addParameter(p, 'WaveDir', waveDir);
    addParameter(p, 'Labels', labelsCsv);
    addParameter(p, 'Only', "");
    parse(p, varargin{:});

    addpath(here);
    arms = {'driven', 'brownian', 'silent'};
    seeds = [101, 202, 303];
    R.waveform = struct();

    for a = 1:numel(arms)
        arm = arms{a};
        family = family_for(arm);
        for s = 1:numel(seeds)
            tag = sprintf('%s_seed%d', arm, seeds(s));
            if strlength(p.Results.Only) > 0 && ~strcmp(tag, char(p.Results.Only))
                continue
            end
            folder = fullfile(p.Results.WaveDir, 'results', ['waveform_' tag]);
            fprintf('\n=== waveform %s ===\n', tag);
            S = classify_waveform(folder, p.Results.Labels, 'Family', family);
            R.waveform.(matlab.lang.makeValidName(tag)) = S;
            if strcmp(arm, 'driven')
                fprintf('--- field-only ---\n');
                Sf = classify_waveform(folder, p.Results.Labels, 'Family', 'field');
                R.waveform.([matlab.lang.makeValidName(tag) '_field']) = Sf;
            end
        end
    end
end

function family = family_for(arm)
    if strcmp(arm, 'brownian')
        family = 'brownian';
    else
        family = 'biology';
    end
end
