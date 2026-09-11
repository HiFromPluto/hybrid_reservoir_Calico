% MATLAB entry point. Run from the HybridDish folder after a BSim job finishes.
%
%   addpath('matlab')
%   S = ridge_narma10('output/driven_seed101');

addpath('matlab');
if ~isfile(fullfile('output', 'driven_seed101', 'voxels.csv'))
    error(['No output/driven_seed101/voxels.csv yet. From this folder run:\n' ...
           '  compile_and_run.cmd config\\driven_seed101.properties']);
end
S = ridge_narma10('output/driven_seed101', 'input/narma10_target.csv');
disp(S);
