%% DFDaddpaths.m

% Script to set paths

path0 = strrep(which('setup.m'),'/setup.m','');

addpath(fullfile(path0,'funcs'));
addpath(fullfile(path0,'aux'));
addpath(genpath(fullfile(path0,'experiments')));
addpath(fullfile(path0,'scripts'));
addpath(fullfile(path0,'figure_scripts'));
clear path0;