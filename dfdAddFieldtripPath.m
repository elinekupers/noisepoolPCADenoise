function dfdAddFieldtripPath(pth)
% Add fieldtrip to Matlab path
% dfdAddFieldtripPath(pth)
%
% Example: 
%   dfdAddFieldtripPath('~/matlab/toolboxes/Fieldtrip/')

    addpath(pth);
    ft_defaults;
return