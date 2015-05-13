% Prepare field trip cfg struct to run by Denoise Wrapper

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

% ******* Load data and design *****************
    tmp = load(sprintf(fullfile(dfdRootPath, 'data', 's0%d_dataAfterICA.mat'),1));
    trial_data = tmp.data.trial;
    
    %%
    n_epochs_per_trial                  = 1;
    nchannels                           = size(trial_data{1}, 1);
    nepochtimepoints                    = 1000;
    ntimepoints_per_trial               = nepochtimepoints * n_epochs_per_trial;
    nepochs                             = n_epochs_per_trial*numel(trial_data); % all runs together

    sensorData = catcell(3, trial_data);
    sensorData = sensorData(:, 1:ntimepoints_per_trial, :);
    sensorData = reshape(sensorData, nchannels, nepochtimepoints, nepochs);

    
    %%
    
    
    
    load(sprintf(fullfile(dfdRootPath, 'data', 's0%d_dataAfterICA_conditions.mat'),1)); 
    tmp = load(sprintf(fullfile(dfdRootPath, 'data', 's0%d_conditions.mat'),1)); conditions = tmp.conditions;

    
    all_trials = find(onsets);
    for ii = 1:size(trl(:,1))
        keepEpochs(ii) = find(all_trials== trl(ii,1));
    
    end
    
    conditions = conditions(keepEpochs);
    
    
    % ******** Make design matrix *****************
    design = zeros(length(conditions), 3);
    design(conditions==1,1) = 1; % condition 1 is full field
    design(conditions==5,2) = 1; % condition 5 is right (??)
    design(conditions==7,3) = 1; % condition 7 is left (??)
    % condition 3 is blank
    
    
    freq = megGetSLandABfrequencies(0:150, 1, 12);

    %  Dummy Denoise for broadband analysis
    optbb.verbose         = 'true';
    optbb.pcchoose        = 0;   
    optbb.preprocessfun   = @hpf;  % preprocess data with a high pass filter for broadband analysis
    evokedfun             = @(x)getstimlocked(x,freq); % function handle to determine noise pool
    evalfun               = @(x)getbroadband(x,freq);  % function handle to compuite broadband

    [results,evalout,denoisedspec,denoisedts] = denoisedata(design,sensorData,evokedfun,evalfun,optbb);
    
    
    fname = sprintf(fullfile(dfdRootPath,'data','s0%d_denoisedData_full'),whichSubject);
    
    parsave([fname '_bb_afterICA.mat'], 'results', results, 'evalout', evalout, ...
        'denoisedspec', denoisedspec, 'denoisedts', denoisedts,...
        'badChannels', badChannels, 'badEpochs', badEpochs, 'opt', optbb)        