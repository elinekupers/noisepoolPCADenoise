function rootPath = DFDrootpath()
% Return the path to the root Denoise Field Data directory
%
% This function must reside in the directory at the base of the
% denoiseFieldData directory structure.  It is used to determine the location of various
% sub-directories.
% 
% Example:
%   fullfile(DFDrootpath,'data')

rootPath=which('DFDrootpath');

rootPath=fileparts(rootPath);

return

%
% This script will set correct paths in MATLAB to use the 'Denoise Field Data'
% algorithm for the paper:
%   AUTHORS. YEAR. TITLE. JOURNAL. VOLUME. ISSUE. DOI.
%
% DFDrootpath.m
%
% Example: 
%   [] = DFDrootpath()


% % Add path with functions
% path0 = strrep(which('DFDrootpath.m'),'/DFDrootpath.m','');
% addpath(genpath(path0));
% 
% addpath(fullfile(path0,'funcs'));
% addpath(fullfile(path0,'aux'));
% addpath(genpath(fullfile(path0,'experiments')));
% addpath(fullfile(path0,'scripts'));
% addpath(fullfile(path0,'figure_scripts'));
% 
% % clear path0;






