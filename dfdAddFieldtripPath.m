function dfdAddFieldtripPath(path)
    if notDefined('path'), path = '/Volumes/server/Projects/MEG/code/fieldtrip'; 
    if ~exist(path,'dir'); path = '/Volumes/server-1/Projects/MEG/code/fieldtrip'; end; end
    addpath(path);
    ft_defaults;
return