function [sensorData, design, badChannels, conditionNames, okEpochs] = ...
    megLoadData(megDataDir,conditionNumbers,opt)
% load meg data from data and format appropriately 
% megDataDir : directory where data are stored (6 matrices plus
%              epoch_conditons)
% conditionNumbers: 1:6 (should be in increasing order) 
%
% TODO : do match conditionNumbers rather than sort to fix ordering
%        requirement
%      : add option to keep bad epochs 
% 

if notDefined('opt'),    opt = struct(); end
if ~isfield(opt,'shuffle_epoch'),        opt.shuffle_epoch = false; end
if ~isfield(opt,'remove_strtend_epoch'), opt.remove_strtend_epoch = false; end
if ~isfield(opt,'badepoch_avgchannum'),  opt.badepoch_avgchannum = 6; end
if ~isfield(opt,'verbose'),              opt.verbose = true; end

if opt.verbose
    fprintf('=============================================================\n');
    fprintf('  Loading data ... \n');
    fprintf('=============================================================\n');
end

% load sensorData
% sensorData = load(fullfile(megDataDir,'data_for_denoising'));
% sensorData = sensorData.data;
conditionNamesAll = {'ON FULL','ON RIGHT','ON LEFT','OFF FULL','OFF RIGHT','OFF LEFT'};
conditionNames    = conditionNamesAll(conditionNumbers);
tepochs    = [];
sensorData = [];
for ii = 1:length(conditionNames)
    dataName = ['ts_', lower(regexprep(conditionNames{ii},' ','_'))];
    %dataName = ['ts_', lower(regexprep(conditionNames{ii},' ','_')), '_epoched'];
    
    if opt.verbose, fprintf('(megLoadData) loading %s\n', dataName); end
    data     = load(fullfile(megDataDir,dataName));
    currdata = data.(dataName);
    tepochs  = cat(1,tepochs,ii*ones(size(currdata,2),1));
    if opt.shuffle_epoch
        currdata = currdata(:,randperm(size(currdata,2)),:);
    end
    sensorData = cat(2,sensorData,currdata);
end

% load conditions
% epoch_conditions:
%   col2: condition numbers as they occurred in the experiment
%   col3: epoch numbers of data matrix as in order in experiment
%   col4: epoch numbers as in order in experiment (added for sanity checking)
epoch_conditions  = load(fullfile(megDataDir,'epoch_conditions'));
epoch_conditions  = epoch_conditions.epochs_condition;
% this gives us a subset of epoch_conditions, just for the conditions we
% want to analyze 
epoch_conditions  = epoch_conditions(ismember(epoch_conditions(:,2),conditionNumbers),:);
epoch_conditions  = [epoch_conditions, (1:size(epoch_conditions,1))'];
% indices for reordering data into experimental order
% now that this require the full dataset loaded to index correctly as the
% indices go up to 1080, and this would be incorrect if we only load a
% subset of the data 
data2condIdx0     = epoch_conditions(:,3); 

% another way of doing the indexing, also used for sanity check
% first, sort by conditions, so now matching data order
[tmp,cond2dataIdx] = sortrows(epoch_conditions,2);
% now sort by epoch order, to get indices for reordering data 
[~, data2condIdx] = sortrows(tmp,4);
% now do our check to see these two things match
% the assertion is only valid for full data set 
if isequal(conditionNumbers,1:6)
    assert(sum(data2condIdx~=data2condIdx0)==0);
end

% arrange data into the conditions as they occurred in the experiment 
sensorData = sensorData(:,data2condIdx,:);

% format data into the right dimensions - chan x time x epochs
sensorData = permute(sensorData,[3,1,2]);
sensorData = sensorData(1:157,:,:);
% replace missing data with nan's
sensorData(sensorData==0) = nan; 

if opt.verbose
    fprintf('number of epochs : %d, number of runs %d\n', size(sensorData,3), size(sensorData,3)/72);
end

okEpochs = true(size(sensorData,3),1);

% remove epochs at the beginning and end of every stimulus presentation 
if opt.remove_strtend_epoch
    okEpochs = reshape(okEpochs,6,[]); % 72 sec runs repeated 15x
    okEpochs([1,6],:) = false;
    okEpochs = okEpochs(:);
end

% remove bad epochs
okEpochs = okEpochs & megIdenitfyBadEpochs(sensorData,0.5);
sensorData = sensorData(:,:,okEpochs);
epoch_conditions = epoch_conditions(okEpochs,:);

% find bad channels - those where at least 50% are nan's 
badChannels = megIdenitfyBadChannels(sensorData, 0.5);
badChannels(98) = 1; % for now we add this in manually
if opt.verbose
    fprintf('\tbadChannels : '); 
    fprintf('%g ', find(badChannels)');
    fprintf('\n');
end

% remove bad channels and replace remaining bad epochs 
net=load('meg160xyz.mat');
sensorData = megReplaceBadEpochs(sensorData,net,[],opt.badepoch_avgchannum);
sensorData = sensorData(~badChannels,:,:);
if opt.verbose, fprintf('\tnumber of final epochs = %d\n', size(sensorData,3)); end

% check that we have no nan's and how many zeros we have
if opt.verbose
    fprintf('(megLoadData) number of data points with 0 : %d\n', sum(vectify(sensorData==0)));
    fprintf('(megLoadData) number of data points with nans : %d\n', sum(vectify(isnan(sensorData))));
    % these are the epochs which averaging neighors still give you nothing
    nanepochs = squeeze(all(sensorData==0,2));
    % nr indexes epochs, nc indexes channels 
    [nr,nc] = find(nanepochs');
    % print the channels and corresponding bad epochs 
    nc2 = unique(nc);
    for jj = 1:length(nc2)
        fprintf('\tchannel %d empty epochs : ', nc2(jj));
        fprintf('%g ', nr(nc==nc2(jj))');
        fprintf('\n');
    end
end

% design matrix
% get condition numbers of all possible ON conditions 
onCondsAll = find(cellfun(@isempty,strfind(conditionNamesAll,'OFF')));
% get the ON conditions for the currently specified conditions 
onConds = conditionNumbers(ismember(conditionNumbers,onCondsAll));
% figure out the dimensions of the design matrix 
design = zeros(size(sensorData,3),length(onConds));
for k = 1:length(onConds)
    design(epoch_conditions(:,2)==onConds(k),k)=1;
end
