%% Doing automatic artifact rejection

% Get a 3D array for the pre, post, relative start point for cfg.trl
project_pth = '/Volumes/server/Projects/MEG/SSMEG';
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
% meg_add_fieldtrip_paths(field_trip_pth, 'yokogawa_defaults');
meg_add_fieldtrip_paths(field_trip_pth);

addpath(genpath('~/matlab/git/meg_utils'));

% Find subjects for this project
subject_pths = dir(fullfile(project_pth,'*SSMEG_*'));

subjects = 1
%% Read in sqd file and preprocess data to get onsets and conditions
data_pth = fullfile(project_pth, subject_pths(subjects).name, 'raw');
[ts, meg_files] = meg_load_sqd_data(data_pth, '*SSMEG_*');

trigger = meg_fix_triggers(ts(:,trigger_channels));
onsets = ssmeg_trigger_2_onsets(trigger, subjects);

[sensorData, conditions] = meg_make_epochs(ts, onsets, epoch_time, fs);

[sensorData, badChannels, badEpochs] = meg_preprocess_data(sensorData(:,:,data_channels), ...
    varThreshold, badChannelThreshold, badEpochThreshold, 'meg160xyz.mat');

%% Define cfg variable for epoching the raw data
cfg = [];

cfg.dataset = fullfile(data_pth,meg_files(subjects).name);
cfg.trialdef.trig = trigger_channels;
cfg.trialdef.eventtype = 'trial';
cfg.trialdef.pre = 0;
cfg.trialdef.post = 1;
cfg.trialdef.trialFunHandle = @ssmeg_ft_trial_fun;
cfg.trialdef.conditions = conditions;
cfg.trialdef.onsets = onsets;

[trl, Events] = ssmeg_ft_trial_fun(cfg, onsets, conditions);
cfg.trl = trl;
cfg.event = Events;

cfg = ft_definetrial(cfg);

%% Preprocess, defining channel as 'MEG', will only take the first 157 channels
%   Takes ~2 minutes
%     cfg.channel = 'MEG';
%     data = ft_preprocessing(cfg);

%     %% Run ICA
%     %   Takes ~20 minutes for whole experiment
%     %   Takes ~ 1 minute for 72 epochs/trials
%     cfg              = [];
%     cfg.method       = 'runica'; % this is the default and uses the implementation from EEGLAB
%     cfg.trials       = 1:72; % all' or a selection given as a 1xN vector (default = 'all')
%     cfg.numcomponent = 'all'; % or number (default = 'all')
%     cfg.demean       = 'no'; % or 'yes' (default = 'yes')
%     comp = ft_componentanalysis(cfg, data);


%% Jump
cfg                    = [];
cfg.trl = trl;
cfg.continuous = 'no';
cfg.datafile = fullfile(data_pth,meg_files(subjects).name);
cfg.headerfile = fullfile(data_pth,meg_files(subjects).name);
data_hdr = load('hdr'); data_hdr = data_hdr.hdr;
cfg.layout = ft_prepare_layout(cfg, data_hdr);

% channel selection, cutoff and padding
cfg.artfctdef.zvalue.channel    = 'MEG';
cfg.artfctdef.zvalue.cutoff     = 20;
cfg.artfctdef.zvalue.trlpadding = 0;
cfg.artfctdef.zvalue.artpadding = 0;
cfg.artfctdef.zvalue.fltpadding = 0;

% algorithmic parameters
cfg.artfctdef.zvalue.cumulative    = 'yes';
cfg.artfctdef.zvalue.medianfilter  = 'yes';
cfg.artfctdef.zvalue.medianfiltord = 9;
cfg.artfctdef.zvalue.absdiff       = 'yes';

% make the process interactive
cfg.artfctdef.zvalue.interactive = 'yes';

[cfg, artifact_jump] = ft_artifact_zvalue(cfg);


% Remove these trials with jump artifact
for ii = 1:size(artifact_jump,1)
    badEpochs1(ii) = find(trl(:,1)==artifact_jump(ii,1));
end

trl(badEpochs1,:,:) = [];



%% Muscle artifact


% muscle
cfg            = [];
cfg.trl        = trl;
cfg.datafile   = fullfile(data_pth,meg_files(subjects).name);
cfg.headerfile = fullfile(data_pth,meg_files(subjects).name);
cfg.continuous = 'no';

% channel selection, cutoff and padding
cfg.artfctdef.zvalue.channel = 'MEG';
cfg.artfctdef.zvalue.cutoff      = 4;
cfg.artfctdef.zvalue.trlpadding  = 0;
cfg.artfctdef.zvalue.fltpadding  = 0;
cfg.artfctdef.zvalue.artpadding  = 0.1;

% algorithmic parameters
cfg.artfctdef.zvalue.bpfilter    = 'yes';
cfg.artfctdef.zvalue.bpfreq      = [110 140];
cfg.artfctdef.zvalue.bpfiltord   = 9;
cfg.artfctdef.zvalue.bpfilttype  = 'but';
cfg.artfctdef.zvalue.hilbert     = 'yes';
cfg.artfctdef.zvalue.boxcar      = 0.2;

% make the process interactive
cfg.artfctdef.zvalue.interactive = 'yes';

[cfg, artifact_muscle] = ft_artifact_zvalue(cfg);


% Remove these trials with muscle and some other weird jumps that got
% trough the first artifact rejection
for ii = 1:size(artifact_muscle,1)
    badEpochs2(ii) = find(trl(:,1)==artifact_muscle(ii,1));
end

trl(badEpochs2,:,:) = [];


%% Slow drift & bad channels
cfg.trl = trl;
cfg.channel = 'MEG';
data = ft_preprocessing(cfg);



cfg          = [];
cfg.method   = 'channel';
cfg.alim     = [];
cfg.megscale = 5;
cfg.keepchannel = 'no';
dummy        = ft_rejectvisual(cfg,data);


%% ICA for eyeblinks, eyemovements, heartbeat.

for ii = 1:length(dummy.label)
    goodChannels(ii) = str2num(dummy.label{ii}(3:end));
end
 % insteaad of 19 or 74)
index69 = find(goodChannels==69)
index99 = find(goodChannels==99)

goodChannels(index69)= []
goodChannels(index99) = []

cfg = [];
cfg.dataset      = fullfile(data_pth,meg_files(subjects).name);
cfg.trl          = trl;
cfg.channel      = goodChannels;
data = ft_preprocessing(cfg);



% Run ICA
%   Takes ~20 minutes for whole experiment
%   Takes ~ 1 minute for 72 epochs/trials
cfg              = [];
cfg.method       = 'runica'; % this is the default and uses the implementation from EEGLAB
cfg.trials       = 'all'; % all' or a selection given as a 1xN vector (default = 'all')
cfg.numcomponent = 'all'; % or number (default = 'all')
cfg.demean       = 'no'; % or 'yes' (default = 'yes')
comp = ft_componentanalysis(cfg, data);

%% Get timelocked averaged events for eyeballing
cfg = [];
timelock = ft_timelockanalysis(cfg, comp);


%% Prepare layout and plot timeseries, components
data_hdr = load('hdr'); data_hdr = data_hdr.hdr;
cfg.layout = ft_prepare_layout(cfg, data_hdr);

% plot the components for visual inspection
figure
cfg = [];
cfg.component = [1:20];       % specify the component(s) that should be plotted
cfg.comment   = 'no';
ft_topoplotIC(cfg, comp)

cfg = [];
cfg.viewmode = 'component';
ft_databrowser(cfg, comp)

%% remove the bad components and backproject the data
cfg = [];

% For Subject 1 remove ICs: [45 51]

cfg.component = [45 51]; % to be removed component(s)
data = ft_rejectcomponent(cfg, comp, data);

fname = fullfile(dfdRootPath, 'data', sprintf('s0%d_dataAfterICA',subjects));
save([fname '.mat'], 'data');

fname = fullfile(dfdRootPath, 'data', sprintf('s0%d_dataAfterICA_conditions',subjects));
save([fname '.mat'], ['trl', 'onsets']);

fname = fullfile(dfdRootPath, 'data', sprintf('s0%d_dataAfterICA_cfg',subjects));
save([fname '.mat'], 'cfg');


