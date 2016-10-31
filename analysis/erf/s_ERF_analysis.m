% s_ERF_analysis()

% Script to plot the eye related function of the data

project_pth                  = '/Volumes/server/Projects/MEG/SSMEG';
subjects                     = 6;
dataChannels                 = 1:157;
erfLength                    = 200; % ms after trigger onset
msConditionOrder             = [3,1,5,7];

%% To run script, you need the Field trip toolbox

% Add fieldtrip path
field_trip_pth = fullfile(fileparts(project_pth), 'code', 'fieldtrip');
% meg_add_fieldtrip_paths(field_trip_pth, 'yokogawa_defaults');
addpath(genpath(field_trip_pth))

% Find subjects for this project
subject_pths = dir(fullfile(project_pth,'*SSMEG_*'));

condNames   = {'Blank','Full','Left','Right'};

dataWithEpochs = {};
dataWithoutEpochs = {};



for whichSubject = subjects
    
    %     ------------------ Load ms timings ----------------------------
    %     tmp = load(sprintf(fullfile(dfdRootPath, 'analysis', 'data', 'eye', 's%02d_epochswithms.mat'),whichSubject)); epochsWithMS = tmp.allEpochsWithMS;
    tmp = load(sprintf(fullfile(dfdRootPath, 'analysis', 'data', 's%02d_conditions.mat'),whichSubject)); conditions = tmp.conditions;
    tmp = load(sprintf(fullfile(dfdRootPath, 'analysis', 'data', 'eye','s0%d_ms.mat'),whichSubject)); ms = tmp.ms;
    tmp = load(sprintf(fullfile(dfdRootPath, 'analysis', 'data', 's%02d_denoisedData_rm1epoch_bb.mat'),whichSubject)); badEpochs = tmp.badEpochs; clear tmp;
    
    nrEpochs = size(cat(1,ms{:}),1);
    % Load sensorData
    tmp = load(sprintf(fullfile(dfdRootPath, 'analysis', 'data', 's%02d_sensorData.mat'),whichSubject)); sensorData = tmp.sensorData; clear tmp;
    sensorData = sensorData(:,1:nrEpochs,:);
    conditions = conditions(1:nrEpochs);
    
    epochLength = size(sensorData,1);
    allDataWithEpochs = {};
    allDataWithoutEpochs = {};
    
    dataWithEpochs = []; dataWithoutEpochs = [];
    
    
    
    % Loop over microsaccade condition order (Blank, Full, Left, Right)
    for cond = msConditionOrder
        
        % Get ms onset times per epoch
        theseMSEpochs = ms{find(cond==msConditionOrder)};
        
        % Get sensordata for this particular condition
        theseData = sensorData(:,conditions==cond,dataChannels);
        
        % Loop over epochs within these MEG and MS data
        for epoch = 1:size(theseMSEpochs,1)
            
            % Loop over epochs for particular conditions
            thesetimePoints = theseData(:,epoch,:);
            
            % Find the timepoint of the MS onset
            msOnsets = find(theseMSEpochs(epoch,:));
            
            % If there is a MS onset
            if ~isempty(msOnsets)
                
                % Loop over onsets in case there are multiple MS in one
                % epoch
                for ii = 1:length(msOnsets);
                    
                    % Find MEG sensorData starting from MS onset + 199
                    % ms.
                    if msOnsets(ii) <= epochLength-erfLength;
                        dataWithEpochs = cat(3,dataWithEpochs, thesetimePoints(msOnsets(ii):msOnsets(ii)+(erfLength-1),:));
                        
                        % If the timepoint is at the end of the epoch (800 ms or later within 1000 ms epoch),
                        % we have to pad with NaNs
                    else
                        indices = thesetimePoints(msOnsets(ii):end,:);
                        indices = [indices;nan(erfLength-size(indices,1),length(dataChannels))];
                        
                        dataWithEpochs = cat(3,dataWithEpochs, indices);
                    end
                end
            else
                % Take a sample of the epoch that doesn't contain a MS
                dataWithoutEpochs = cat(3,dataWithoutEpochs, thesetimePoints(1:200,:));
            end
            
        end
        
        % dataWithEpochs or dataWithoutEpochs has dimensions: Timepoints x channels x epoch
        allDataWithEpochs{cond} = dataWithEpochs; dataWithEpochs = [];
        allDataWithoutEpochs{cond} = dataWithoutEpochs; dataWithoutEpochs = [];
        
    end
    
    allDataWithEpochs = allDataWithEpochs(~cellfun('isempty',allDataWithEpochs));
    allDataWithoutEpochs = allDataWithoutEpochs(~cellfun('isempty',allDataWithoutEpochs));
    
    
    %% Visualize ERFs
    cmap = parula(4);
    for chan = dataChannels
        
        figure(1); clf; set(gcf,'Name',sprintf('Channel %d', chan), 'Color','w')
        hold all;
        
        for ii = 1:length(msConditionOrder)
            
            subplot(4,1,ii)
            plot(0:199, nanmean(allDataWithEpochs{ii}(:,chan,:),3)', 'Color', cmap(ii,:), 'LineWidth',2); hold on;
            plot(0:199, nanmean(allDataWithoutEpochs{ii}(:,chan,:),3)', 'Color', cmap(ii,:), 'LineStyle','--', 'LineWidth',2);
            plot(0:199, zeros(1,200), 'k');
            title(condNames(ii));
            legend('With MS', 'Without MS')
            box off
            ylim([-300,300]);
            
        end
%         pause(1)
        hgexport(1,sprintf(fullfile(dfdRootPath, 'analysis', 'figures','erf', 's%02d_chan%d.mat'),whichSubject, chan));
    end
    
end

