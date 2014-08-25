% BULK PREPROCESS ALL SESSION DATA
% SEE ALSO denoisescript_meg.m
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Session Numbers - see megGetDataPaths for what these correspond to
% Note these numbers are just in the order that Helena processed the data,
% and doesn't correspond to any logical way of organizing the sessions.
% sessions 1 and 2 are replaced by sessions 11 and 12, respectively.
% sessions 7 and 8 correspond to attention sessions and should not analyzed
% together with the other sessions
% for manuscript purposes, sessions analyzed should be [11:12,3:6,9:10]
% In the future, session numbers can be re-arranged into something
% simpler (like 1:8), as long as correspondence between old numbers and
% new numbers can be kept track of nicely.
sessionNums = 1:12;

% In mm. Corresponds to a parameter in megLoadData, or the neighborhood for
% doing weighted averaging so that missing data for a given channel can be
% replaced by the weighted average of its neighbors.
% This corresponds to the values Helena used for each of the sessions for
% majority of the analyses.
% TODO: this should be replaced by a weighted average of all other
% channels, rather than those within a neighborhood
badEpochAvg = [6,6,6,8,8,6,8,8,6,6,8,8];

% A string attached to the end of session name to distinguish the kind of
% processing that's been done. see Readme.txt in inputDataDir
sensorDataStr = 'b2';

% Whether to save data
saveData = false;
% Whether to save epoch group data
saveEpochGroup = false;

% data directory to save things
inputDataDir = '/Volumes/HelenaBackup/denoisesuite/tmpmeg/';
% condition Numbers
conditionNumbers = 1:6;

% process the sessions we request 
for ii = sessionNums
    % get session name and top directory
    % this can be modified so that top directory points somewhere else
    [dataset,megDataDir] = megGetDataPaths(ii,conditionNumbers);
    % or we can hard code it like so
    % megDataDir = '/Volumes/server/Projects/MEG/SSMEG/';
    
    % get directory for current session
    megDataDir = fullfile(megDataDir,dataset);
    disp(dataset);
    
    % load session data
    clear loadopt
    loadopt.badepoch_avgchannum  = badEpochAvg(ii);
    [sensorData, design, badChannels, conditionNames, okEpochs] ...
        = megLoadData(megDataDir,conditionNumbers,loadopt);
    
    % save preprocessed data, if requested
    if saveData
        save(fullfile(inputDataDir,'inputdata',sprintf('%s%s',dataset,sensorDataStr)),'sensorData', 'design', 'badChannels');
        fprintf('saved : %s\n', dataset);
    end
    
    % create epoch groups and save them, if requested
    group_epoch = 6;
    shift_epochs = [0,3];
    remove_unequal_group = true;
    % create two, one with shifted epoch groups, one without
    for jj = 1:2
        epochGroup = megEpochGroup(okEpochs,group_epoch,shift_epochs(jj),remove_unequal_group);
        assert(length(epochGroup)==size(sensorData,3)); % sanity check
        fname = sprintf('%s%s_epochGroup%d',dataset,sensorDataStr,group_epoch);
        if shift_epochs(jj), fname = [fname,'s']; end
        if remove_unequal_group, fname = [fname,'o']; end
        if saveEpochGroup
            disp(fname);
            save(fullfile(inputDataDir, 'epochGroups', fname),'epochGroup');
        end
    end
    
end