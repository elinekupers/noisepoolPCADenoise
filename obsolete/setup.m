path0 = strrep(which('setup.m'),'/setup.m','');

% EXTERNAL TOOLBOXES: please change paths accordiningly 
% ------------------------------------------------------

% external toolboxes required for meg 
warning on;
path00 = fileparts(path0);
megtool = {'fieldtrip'};
for k = 1:length(megtool)
    megtoolpath = fullfile(path00,megtool{k}); %<-- UPDATE HERE
    if ~exist(megtoolpath,'dir')
        warning('denoisesuite:setup','toolbox %s not found! required for MEG %s', megtool{k},megtoolpath),
        megtoolpath = uigetdir(path00);
    end
    if strcmpi('fieldtrip', megtoolpath)
        addpath(megtoolpath)
        run('ft_defaults');
    else
        addpath(genpath(megtoolpath));
    end
end
    
% external toolboxes required for eeg
eegtool = {'CSDtoolbox','svndl_code'};
for k = 1:length(eegtool)
    eegtoolpath = fullfile(path00,'EEG',eegtool{k}); %<-- UPDATE HERE
    if ~exist(eegtoolpath,'dir')
        warning('denoisesuite:setup','toolbox %s not found! required for EEG %s', eegtool{k},eegtoolpath),
        eegtoolpath = uigetdir(path00);
    else
        addpath(genpath(eegtoolpath));
        eegtoolpath = uigetdir(path00);
    end
end


addpath(fullfile(path0,'funcs'));
addpath(fullfile(path0,'aux'));
addpath(genpath(fullfile(path0,'experiments')));
addpath(fullfile(path0,'scripts'));
addpath(fullfile(path0,'figure_scripts'));
clear;