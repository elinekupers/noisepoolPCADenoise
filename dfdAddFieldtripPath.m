function dfdAddFieldtripPath(path)
    if notDefined('path'), path = '/Volumes/server/Projects/MEG/code/fieldtrip'; end
    addpath(path);
    ft_defaults;
return