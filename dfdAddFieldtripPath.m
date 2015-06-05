function dfdAddFieldtripPath(path)
    if notDefined('path'), path = '/Volumes/server/Projects/MEG/code/fieldtrip'; end
    addpath(genpath(path));
    ft_defaults;
return