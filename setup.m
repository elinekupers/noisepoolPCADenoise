path0 = strrep(which('setup.m'),'/setup.m','');
addpath(fullfile(path0,'funcs'));
addpath(fullfile(path0,'aux'));

% this is temporary - these functions will be moved 
path00 = fileparts(path0);
addpath(fullfile(path00,'denoisecode'))
addpath(genpath(fullfile(path00,'fieldtrip')))
