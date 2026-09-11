function [X, names, u, windows] = window_features(T, family)
%WINDOW_FEATURES  Collapse 16 intra-window samples to one row per window.
%   State channels (R, L, AHL, Den) are window means.
%   Count channels (Input_Driven_Death) use the last sample of the window.
%
%   family:
%     'biology'  — Receiver_R + Lum_Mean + Input_Driven_Death (408 cols)
%     'lum'      — Lum_Mean_* window mean (200 cols; KR diagnostic)
%     'receiver' — Receiver_R_* window mean (200 cols; KR diagnostic)
%     'brownian' — Den_* (200 cols)
%     'field'    — AHL_uM_* (200 cols; diagnostic, not the official readout)

    if nargin < 2 || isempty(family)
        family = 'biology';
    end
    if ~ismember('Window', T.Properties.VariableNames)
        error('window_features:header', 'voxels table needs a Window column');
    end

    switch lower(family)
        case 'biology'
            meanPref = {'Receiver_R_', 'Lum_Mean_'};
            lastPref = {'Input_Driven_Death_'};
        case {'lum', 'lum_mean'}
            meanPref = {'Lum_Mean_'};
            lastPref = {};
        case {'receiver', 'receiver_r'}
            meanPref = {'Receiver_R_'};
            lastPref = {};
        case 'brownian'
            meanPref = {'Den_'};
            lastPref = {};
        case 'field'
            meanPref = {'AHL_uM_'};
            lastPref = {};
        otherwise
            error('window_features:family', 'Unknown family: %s', family);
    end

    allNames = T.Properties.VariableNames;
    meanNames = names_with_prefix(allNames, meanPref);
    lastNames = names_with_prefix(allNames, lastPref);
    names = [meanNames, lastNames];
    if isempty(names)
        error('window_features:empty', 'No %s columns in this voxels table', family);
    end

    windows = unique(T.Window, 'stable');
    nW = numel(windows);
    X = zeros(nW, numel(names));
    u = zeros(nW, 1);

    hasU = ismember('Input_AC1_AHL', allNames);
    for i = 1:nW
        rows = T.Window == windows(i);
        block = T(rows, :);
        if ~isempty(meanNames)
            X(i, 1:numel(meanNames)) = mean(block{:, meanNames}, 1, 'omitnan');
        end
        if ~isempty(lastNames)
            X(i, numel(meanNames)+1:end) = block{end, lastNames};
        end
        if hasU
            u(i) = block.Input_AC1_AHL(end);
        end
    end
end

function out = names_with_prefix(allNames, prefixes)
    out = {};
    for p = 1:numel(prefixes)
        hit = allNames(startsWith(allNames, prefixes{p}));
        out = [out, hit]; %#ok<AGROW>
    end
end
