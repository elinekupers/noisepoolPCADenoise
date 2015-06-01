%% Script to run ICA on SSMEG data, with field trip toolbox, visualize
% the different components and possibility to reject 


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

for subjects = 1
    %% Read in sqd file and preprocess data to get onsets and conditions
    data_pth = fullfile(project_pth, subject_pths(subjects).name, 'raw');
    [ts, meg_files] = meg_load_sqd_data(data_pth, '*SSMEG_*');
    
    trigger = meg_fix_triggers(ts(:,trigger_channels));
    onsets = ssmeg_trigger_2_onsets(trigger, subjects);
    
    [sensorData, conditions] = meg_make_epochs(ts, onsets, epoch_time, fs);
    
%     [sensorData, badChannels, badEpochs] = meg_preprocess_data(sensorData(:,:,data_channels), ...
%         varThreshold, badChannelThreshold, badEpochThreshold, 'meg160xyz.mat');
    
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
    cfg.channel = 'MEG';
    data = ft_preprocessing(cfg);
    
    %% Run ICA
    %   Takes ~20 minutes for whole experiment
    %   Takes ~ 1 minute for 72 epochs/trials
    cfg              = [];
    cfg.method       = 'runica'; % this is the default and uses the implementation from EEGLAB
    cfg.trials       = 1:72; % all' or a selection given as a 1xN vector (default = 'all')
    cfg.numcomponent = 'all'; % or number (default = 'all')
    cfg.demean       = 'no'; % or 'yes' (default = 'yes')
    comp = ft_componentanalysis(cfg, data);
    
    %% Get timelocked averaged events for eyeballing
    cfg = [];
    timelock = ft_timelockanalysis(cfg, comp);

    % Or either get the Fourier spectrum
%     cfg = [];
%     cfg.method = 'wavelet';
%     freq = ft_freqanalysis(cfg, comp);
    
    %% Prepare layout and plot timeseries, components
    data_hdr = load('hdr'); data_hdr = data_hdr.hdr;
    cfg.layout = ft_prepare_layout(cfg, data_hdr);

    % plot the components for visual inspection
    figure
    cfg = [];
    cfg.component = [21:40];       % specify the component(s) that should be plotted
    cfg.comment   = 'no';
    ft_topoplotIC(cfg, comp)

    cfg = [];
    cfg.viewmode = 'component';
    ft_databrowser(cfg, comp)

    %% remove the bad components and backproject the data
    cfg = [];
    
    % For Subject 1 remove ICs: [11 16 22]
    
    cfg.component = []; % to be removed component(s) 
    data = ft_rejectcomponent(cfg, comp, data);
    
    fname = fullfile(dfdRootPath, 'data' sprintf('s0%d_dataAfterICA',subjects));
    save([fname '.mat'], 'data');
end