function [sensorData, design, badChannels, epochGroup] = ...
    megLoadData(megDataDir,conditionNumbers,opt)
% load meg data from data and format appropriately 
% conditionNumbers: 1 FULL, 2 RIGHT, 3 LEFT, 0 BLANK

if notDefined('opt'),    opt = struct(); end
if ~isfield(opt,'group_epoch'),   opt.group_epoch  = 1;  end
%if ~isfield(opt,'shuffle_epoch'), opt.shuffle_epoch = false; end
if ~isfield(opt,'remove_strtend_epoch'), opt.remove_strtend_epoch = false; end
if ~isfield(opt,'badepoch_avgchannum'),  opt.badepoch_avgchannum = 3; end
if ~isfield(opt,'verbose'),       opt.verbose = true; end

% load sensorData and conditions 
sensorData = load(fullfile(megDataDir,'data_for_denoising'));
sensorData = sensorData.data;
epoch_conditions = load(fullfile(megDataDir,'epoch_conditions'));
epoch_conditions = epoch_conditions.epochs_condition_count;
idx = ismember(epoch_conditions(:,2),conditionNumbers); 

% analyze only specified conditions 
sensorData = sensorData(:,idx,:);
epoch_conditions = epoch_conditions(idx,:);

% format data into the right dimensions - chan x time x epochs
sensorData = permute(sensorData,[3,1,2]);
sensorData = sensorData(1:157,:,:);
% replace missing data with nan's
sensorData(sensorData==0) = nan; 

% group by trial
if opt.group_epoch > 1
    if opt.verbose
        fprintf('(megLoadData) grouping everying %d epochs for denoising\n', opt.group_epoch);
    end
    epochGroup = repmat(1:size(sensorData,3)/opt.group_epoch,opt.group_epoch,1);
    epochGroup = epochGroup(:);
else
    epochGroup = [];
end

if opt.verbose
    fprintf('number of epochs : %d, number of runs %d\n', size(sensorData,3), size(sensorData,3)/72);
end

% remove epochs at the beginning and end of every stimulus presentation 
if opt.remove_strtend_epoch
    kept_epochs = true(size(sensorData,3),1);
    kept_epochs = reshape(kept_epochs,6,[]); % 72 sec runs repeated 15x
    kept_epochs([1,6],:) = false;
    kept_epochs = kept_epochs(:);
    sensorData = sensorData(:,:,kept_epochs);
    epoch_conditions = epoch_conditions(:,:,kept_epochs);
    epochGroup = epochGroup(kept_epochs);
end

% remove bad epochs
[sensorData,okEpochs] = megRemoveBadEpochs({sensorData},0.5);
epoch_conditions = epoch_conditions(okEpochs{1},:);
if ~notDefined('epochGroup'), epochGroup = epochGroup(okEpochs{1}); end

% find bad channels - those where at least 50% are nan's 
badChannels = megIdenitfyBadChannels(sensorData, 0.5);
badChannels(98) = 1; % for now we add this in manually
if opt.verbose, fprintf('\tbadChannels = %g \n', find(badChannels)'); end

% remove bad epochs and channels
net=load('meg160xyz.mat');
sensorData = megReplaceBadEpochs(sensorData{1},net,[],opt.badepoch_avgchannum);
sensorData = sensorData(~badChannels,:,:);
if opt.verbose, fprintf('\tnumber of epochs = %d\n', size(sensorData,3)); end

% check that we have no nan's and how many zeros we have
if opt.verbose
    fprintf('(megLoadData) number of chan with 0 : %d\n', sum(vectify(sensorData==0)));
    fprintf('(megLoadData) number of chan with nans : %d\n', sum(vectify(isnan(sensorData))));
    % these are the epochs which averaging neighors still give you nothing
    nanepochs = squeeze(all(sensorData==0,2));
    [nr,nc] = find(nanepochs');
    nc2 = unique(nc);
    for jj = 1:length(nc2)
        fprintf('channel %d epochs:\n', nc2(jj));
        disp(nr(nc==nc2(jj))');
        fprintf('\n');
    end
end

% design matrix
onConds = conditionNumbers(conditionNumbers~=0);
design = zeros(size(sensorData,3),nnz(conditionNumbers));
for k = 1:length(onConds)
    design(epoch_conditions(:,2)==onConds(k),k)=1;
end