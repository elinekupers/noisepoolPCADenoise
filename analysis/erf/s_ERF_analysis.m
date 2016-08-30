% s_ERF_analysis()

% Script to plot the eye related function of the data

project_pth                   = '/Volumes/server/Projects/MEG/SSMEG';
subjects                      = 6:8;
dataChannels                 = 1:157;
triggerChannels              = 161:164;
fs                            = 1000;        % sample rate (Hz)
epoch_time                    = [0 1];       % start and end of epoch, relative to onset, in seconds
%% To run script, you need the Field trip toolbox

% Add fieldtrip path
field_trip_pth = fullfile(fileparts(project_pth), 'code', 'fieldtrip');
% meg_add_fieldtrip_paths(field_trip_pth, 'yokogawa_defaults');
addpath(genpath(field_trip_pth))

% Find subjects for this project
subject_pths = dir(fullfile(project_pth,'*SSMEG_*'));

condNames   = {'Blank','Full','Left','Right'};
epochsInExpWithMS = [];

for whichSubject = subjects
    
    %     ------------------ Load ms timings ----------------------------
    tmp = load(sprintf(fullfile(dfdRootPath, 'analysis', 'data', 'eye', 's%02d_epochswithms.mat'),whichSubject)); epochsWithMS = tmp.allEpochsWithMS;
    tmp = load(sprintf(fullfile(dfdRootPath, 'analysis', 'data', 's%02d_conditions.mat'),whichSubject)); conditions = tmp.conditions;
    
    triggerNrs = unique(tmp.conditions);
    triggerNrs = [triggerNrs(2),triggerNrs(1),triggerNrs(3),triggerNrs(4)]; % 3 = blank, 1 = full, 5 = left, 7 = right
    
    conditionsOrdered = [conditions, [1:length(conditions)]'];
    [cond, idx] = sort(conditionsOrdered,1);
    
    for ii = 1:length(triggerNrs)
        
        epochNrInOrder = find(cond(:,1)==triggerNrs(ii));
        epochNrInExperiment = idx(epochNrInOrder,1);
        epochsInExpWithMS = [epochsInExpWithMS; epochNrInExperiment(epochsWithMS{ii})];
        
    end
    
    save(sprintf(fullfile(dfdRootPath, 'analysis', 'data', 's%02d_msepochidx.mat'),whichSubject), 'epochsInExpWithMS');

    
    % Load sensorData
%     tmp = load(sprintf(fullfile(dfdRootPath, 'analysis', 'data', 's%02d_sensorData.mat'),whichSubject)); sensorData = tmp.sensorData;
    
    %
end

