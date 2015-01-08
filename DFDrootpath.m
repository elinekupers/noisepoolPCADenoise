function path0 = DFDrootpath()
%
% This script will set correct paths in MATLAB to use the 'Denoise Field Data'
% algorithm for the paper:
%   AUTHORS. YEAR. TITLE. JOURNAL. VOLUME. ISSUE. DOI.
%
% DFDrootpath.m
%
% Example: 
%   [] = DFDrootpath()


% Add path with functions
path0 = strrep(which('DFDrootpath.m'),'/DFDrootpath.m','');
addpath(genpath(path0));

addpath(fullfile(path0,'funcs'));
addpath(fullfile(path0,'aux'));
addpath(genpath(fullfile(path0,'experiments')));
addpath(fullfile(path0,'scripts'));
addpath(fullfile(path0,'figure_scripts'));

% clear path0;






