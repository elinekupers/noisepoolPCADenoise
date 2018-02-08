function rootPath = nppRootPath()
% Return the path to the root Noisepool-PCA algorithm directory
%
% This function must reside in the directory at the base of the
% Noisepool-PCA directory structure.  It is used to determine the location of various
% sub-directories.
% 
% Example:
%   fullfile(nppRootpath,'data')

rootPath=which('nppRootPath');

rootPath=fileparts(rootPath);

return
