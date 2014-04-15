path0 = strrep(which('setup.m'),'/setup.m','');
addpath(fullfile(path0,'funcs'));
addpath(fullfile(path0,'aux'));

% this is temporary - these functions will be moved so the toolbox is
% standalone
path00 = fileparts(path0);
addpath(fullfile(path00,'denoisecode'))
addpath(genpath(fullfile(path00,'fieldtrip')))

addpath('/Users/Helena/Work/EEG/EEG/Code');
addpath(genpath('/Users/Helena/Work/EEG/EEG/CSDtoolbox'));
addpath(genpath('/Users/Helena/Work/EEG/EEG/svndl_code'));
