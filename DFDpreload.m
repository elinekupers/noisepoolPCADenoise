function [sensorData, design, badChannels, conditionNames, okEpochs] = ...
                DFDpreload(sessionNums, sensorDataStr, saveData, saveEpochGroup, ...
                            inputDataDir, conditionNumbers)
%% Loads and preproccesses sample MEG data sets before it will be denoised 
%  by the 'Denoise Field Data' algorithm for the paper:
%
%   AUTHORS. YEAR. TITLE. JOURNAL. VOLUME. ISSUE. DOI.
%
% [sensorData, design, badChannels, conditionNames, okEpochs] = ...
%                DFDpreload(sessionNums, sensorDataStr, saveData, saveEpochGroup, ...
%                           inputDataDir, conditionNumbers)
%
% Inputs
%   sessionNums:        Vector of data sets (1 to 8)
%                       [default=1:8]
%   sensorDataStr:      A string attached to the end of session name to distinguish the kind of
%                       processing that's been done. see Readme.txt in inputDataDir
%                       [default='b2']
%   saveData:           Boolean whether to save preprocessed datasets or not
%                       [default=true]
%   saveEpochGroup:     Boolean Whether to save shuffled epoch group data (this 
%                       will probably be taken out of this function)
%                       [default=false]
%   inputDataDir:       String defining data directory to save data
%                       [default=fullfile(DFDrootpath, 'data')]
%   conditionNumbers:   Vector of conditions to use (1 to 6)
%                       [default=1:6]
%
% Outputs
%   sensorData:         Matrix containing data (channels x timepoints x nr epochs)
%   design:             Matrix containing condition design (nr epochs x nr conditions)
%   badChannels:        Vector whether channels are defined as 'bad' (1=good,0=bad)
%   conditionNames:     Vector of conditions used (1 to 6)
%   okEpochs:           Vector whether channels are defined as 'ok' (1=good,0=bad)
%
% Example: Preload 1 subject condition full on and off, save data
%   [sensorData, design, badChannels, conditionNames, okEpochs] = ...
%                DFDpreload(1, [], [], [], [], [1,4])


%% Comments from Helena
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
% sessionNums = 1:12;

% In mm. Corresponds to a parameter in megLoadData, or the neighborhood for
% doing weighted averaging so that missing data for a given channel can be
% replaced by the weighted average of its neighbors.
% This corresponds to the values Helena used for each of the sessions for
% majority of the analyses.

% TODO: this should be replaced by a weighted average of all other
% channels, rather than those within a neighborhood
% badEpochAvg = [6,6,6,8,8,6,8,8,6,6,8,8]; % Helena's order
badEpochAvg = [6,8,8,6,6,6,8,8]; % Our datasets order

% A string attached to the end of session name to distinguish the kind of
% processing that's been done. see Readme.txt in inputDataDir
% sensorDataStr = 'b2';

% Whether to save data
% saveData = true;
% Whether to save epoch group data
% saveEpochGroup = true;

% data directory to save things
% inputDataDir = '/Volumes/server/Projects/MEG/GLMdenoised/tmpmeg';
% condition Numbers
% conditionNumbers = 1:6;
%% Define options in case they are not specified by the user

if notDefined('opt'),    opt = struct(); end
if ~isfield(opt,'sessionNums'),          opt.sessionNums = 1:8; end % Do all subjects
if ~isfield(opt,'sensorDataStr'),        opt.sensorDataStr = 'b2'; end % 
if ~isfield(opt,'saveData'),             opt.saveData = true; end %
if ~isfield(opt,'saveEpochGroup'),       opt.saveEpochGroup = true; end %
if ~isfield(opt,'inputDataDir'),         opt.inputDataDir = fullfile(DFDrootpath, 'data'); end %
if ~isfield(opt,'conditionNumbers'),     opt.conditionNumbers = 1:6; end %

%% Process the sessions we request 
for ii = sessionNums
    % Get session name and top directory
    % This can be modified so that top directory points somewhere else
    [dataset,DFDDataDir] = DFDgetdatapaths(ii,conditionNumbers,inputDataDir);
    
    % Get directory for current session, or point to download function
    try
        DFDDataDir = fullfile(DFDDataDir,dataset);
    catch ME
        idSegLast = regexp(ME.identifier, 'Undefined');
        if idSegLast
            fprintf('(DFDPreLoad): ERROR. Data cannot be found. Please download data with the function DFDdownloaddata');
        else 
            disp(ME)
        end
    end
    
    disp(dataset);
    
    % load session data
    clear loadopt
    loadopt.badepoch_avgchannum  = badEpochAvg(ii);
    [sensorData, design, badChannels, conditionNames, okEpochs] ...
        = megLoadData(DFDDataDir,conditionNumbers,loadopt);
    
    % save preprocessed data, if requested
    if saveData
        if ~exist(fullfile(inputDataDir,'savedProcData'),'dir'); mkdir(fullfile(inputDataDir,'savedProcData')); end
        save(fullfile(inputDataDir,'savedProcData',sprintf('%s%s',dataset,sensorDataStr)),'sensorData', 'design', 'badChannels','okEpochs');
        fprintf('saved : %s\n', dataset);
    end
    
    %% TODO Think about getting this out of this function
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