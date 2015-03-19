function rootPath = dfdRootPath()
% Return the path to the root Denoise Field Data directory
%
% This function must reside in the directory at the base of the
% denoiseFieldData directory structure.  It is used to determine the location of various
% sub-directories.
% 
% Example:
%   fullfile(DFDrootpath,'data')

rootPath=which('dfdRootPath');

rootPath=fileparts(rootPath);

return
