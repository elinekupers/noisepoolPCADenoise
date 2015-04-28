%% Prepare data to perform ICA on example subject

% Run s_SSMEG_prepare_data only for subject 1
% load denoised data set for subject 1 to get Bad


% Define variables
project_pth                   = '/Volumes/server/Projects/MEG/SSMEG';
subjects                      = 1; 
data_channels                 = 1:157;
trigger_channels              = 161:164;
photodiode_channel            = 192; 

fs                            = 1000;        % sample rate (Hz)
epoch_time                    = [0 1];       % start and end of epoch, relative to onset, in secondsfield_trip_pth

varThreshold                  = [0.05 20];  
badChannelThreshold           = 0.2;       
badEpochThreshold             = 0.2;

% Add fieldtrip path
field_trip_pth = fullfile(fileparts(project_pth), 'code', 'fieldtrip');
meg_add_fieldtrip_paths(field_trip_pth, 'yokogawa_defaults');

% Add meg_utils path after field trip path in order to overwrite the
% mytrialfun function. (In case that doesn't work, replace fieldtrip
% function name with an extra character at the end, so there will be only
% one version of this function.
addpath(genpath('~/matlab/git/meg_utils'));

% IMPORTANT:
% Check: which mytrialfun_all.m  is equal to /Users/winawerlab/matlab/git/meg_utils/utils/mytrialfun_all.m

% Change the private field trip function "read_trigger.m" 
% Which is in "/Volumes/server/Projects/MEG/code/fieldtrip/fileio/private/read_trigger.m"  % Private to fileio
% into:
% "/Volumes/server/Projects/MEG/code/fieldtrip/fileio/private/read_trigger0.m",
% so adding a zero at the end of the filename. That is because the private
% functions cannot not be replaced by adding a function with addpath().
% So next check is: which read_trigger.m is equal to Users/winawerlab/matlab/git/meg_utils/ft_development/read_trigger.m


% Find subjects for this project
subject_pths = dir(fullfile(project_pth,'*SSMEG_*'));


% In order to get badChannels, read in sqd file and preprocess data
data_pth = fullfile(project_pth, subject_pths(subjects).name, 'raw');
    [ts, meg_files] = meg_load_sqd_data(data_pth, '*SSMEG_*');
    
trigger = meg_fix_triggers(ts(:,trigger_channels));
onsets = ssmeg_trigger_2_onsets(trigger, subjects);
    
[sensorData, conditions] = meg_make_epochs(ts, onsets, epoch_time, fs);

[sensorData, badChannels, badEpochs] = meg_preprocess_data(sensorData(:,:,data_channels), ...
    varThreshold, badChannelThreshold, badEpochThreshold, 'meg160xyz.mat');

% Define cfg variable
cfg = [];

cfg.dataset = fullfile(data_pth,meg_files(subjects).name);
cfg.badChannels = find(badChannels);
cfg.trialdef.nTrigsExpected = 1080;
cfg.trialdef.trig = trigger_channels;
cfg.trialdef.prestim = 0;
cfg.trialdef.poststim = 1;
% cfg.trialdef.trialFunHandle = @mytrialfun_all(cfg, 50, cfg.trialdef.nTrigsExpected);
cfg.trialdef.trialFunHandle = @mytrialfun_all;



% Use meg utils PCA ICA function
[ft_cleandata, WorkFlow] = meg_pca_ica(cfg.dataset, badChannels, cfg.trialdef);