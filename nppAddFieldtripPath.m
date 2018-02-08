function nppAddFieldtripPath(pth)
% Add fieldtrip to Matlab path
% nppAddFieldtripPath(pth)
%
% Example: 
%   nppAddFieldtripPath('~/matlab/toolboxes/Fieldtrip/')

    addpath(pth);
    ft_defaults;
return