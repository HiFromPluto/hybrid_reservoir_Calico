function S = ridge_narma10(voxel_dir, target_csv, varargin)
%RIDGE_NARMA10  Official 408-feature NARMA-10 readout for HybridDish.
%
%   S = RIDGE_NARMA10('output/driven_seed101')
%   S = RIDGE_NARMA10(voxel_dir, 'input/narma10_target.csv')
%   S = RIDGE_NARMA10(..., 'Family', 'biology')  % or 'brownian' / 'field'
%
%   Washout 40 / train 110 / test 50. Lambda is fit on windows 128..149 only.

    if nargin < 2 || isempty(target_csv)
        here = fileparts(mfilename('fullpath'));
        target_csv = fullfile(fileparts(here), 'input', 'narma10_target.csv');
    end

    p = inputParser;
    addParameter(p, 'Family', 'biology');
    parse(p, varargin{:});

    T = load_voxels(voxel_dir);
    [X, names, u_vox] = window_features(T, p.Results.Family);
    target = readtable(target_csv, 'Delimiter', ';', 'FileType', 'text');
    y = target.y_next;
    u = target.u;

    if numel(y) ~= size(X, 1)
        error('ridge_narma10:windows', ...
            'Target has %d rows, voxels have %d windows', numel(y), size(X, 1));
    end
    if ~isempty(u_vox) && max(abs(u_vox(:) - u(:))) > 1e-6
        warning('ridge_narma10:u', 'Voxel Input_AC1_AHL does not match target u');
    end

    S = ridge_closed(X, y);
    S.family = p.Results.Family;
    S.n_features = numel(names);
    S.feature_names = names;
    S.voxel_dir = voxel_dir;

    fprintf('HybridDish NARMA-10 (%s, %d features)\n', S.family, S.n_features);
    fprintf('  lambda      = %.4g\n', S.lambda);
    fprintf('  train NRMSE = %.4f\n', S.train_nrmse);
    fprintf('  test  NRMSE = %.4f\n', S.test_nrmse);
    fprintf('  test  R^2   = %.4f\n', S.test_r2);
end
