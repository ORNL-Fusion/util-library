function setup_matlab_paths(util_root)
% setup_matlab_paths(util_root)
%
% Add util-library MATLAB routines to the active MATLAB path.

if nargin < 1 || isempty(util_root)
    util_root = fileparts(mfilename('fullpath'));
end

addpath(genpath(fullfile(util_root,'matlab','bfield_library_jdl')))
fprintf('util-library MATLAB paths enabled: %s\n', util_root)
