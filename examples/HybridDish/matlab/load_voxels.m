function T = load_voxels(path)
%LOAD_VOXELS  Read HybridDish voxels.csv (semicolon-delimited, like reservoir_new).
%   T = LOAD_VOXELS(path) accepts a folder that contains voxels.csv, or the
%   CSV path itself. Columns keep their original names.

    if isfolder(path)
        path = fullfile(path, 'voxels.csv');
    end
    if ~isfile(path)
        error('load_voxels:missing', 'No voxels.csv at %s', path);
    end

    opts = detectImportOptions(path, 'Delimiter', ';', 'FileType', 'text');
    opts = setvartype(opts, 'double');
    T = readtable(path, opts);
    T.Properties.VariableNames = matlab.lang.makeValidName( ...
        strtrim(T.Properties.VariableNames));
end
