function [sensorData, design, badChannels, conditionNames, okEpochs] = ...
    megLoadData(megDataDir,conditionNumbers,opt)
% load meg data from data and format appropriately 

if notDefined('opt'),    opt = struct(); end
if ~isfield(opt,'shuffle_epoch'),        opt.shuffle_epoch = false; end
if ~isfield(opt,'remove_strtend_epoch'), opt.remove_strtend_epoch = false; end
if ~isfield(opt,'badepoch_avgchannum'),  opt.badepoch_avgchannum = 3; end
if ~isfield(opt,'verbose'),              opt.verbose = true; end

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
epoch_conditions = load(fullfile(megDataDir,'epoch_conditions'));
epoch_conditions = epoch_conditions.epochs_condition;
epoch_conditions = [epoch_conditions, (1:size(epoch_conditions,1))'];
epoch_conditions = epoch_conditions(ismember(epoch_conditions(:,2),conditionNumbers),:);
data2condIdx = epoch_conditions(:,3);

% sanity check
[tmp,cond2dataIdx2] = sortrows(epoch_conditions,2);
[~,  data2condIdx2] = sortrows(tmp,4);
assert(sum(data2condIdx2~=data2condIdx)==0);

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
if opt.verbose, fprintf('\tbadChannels = %g \n', find(badChannels)'); end

% remove bad channels and replace remaining bad epochs 
net=load('meg160xyz.mat');
sensorData = megReplaceBadEpochs(sensorData,net,[],opt.badepoch_avgchannum);
sensorData = sensorData(~badChannels,:,:);
if opt.verbose, fprintf('\tnumber of epochs = %d\n', size(sensorData,3)); end

% check that we have no nan's and how many zeros we have
if opt.verbose
    fprintf('(megLoadData) number of data points with 0 : %d\n', sum(vectify(sensorData==0)));
    fprintf('(megLoadData) number of data points with nans : %d\n', sum(vectify(isnan(sensorData))));
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
onConds = find(cellfun(@isempty,strfind(conditionNames,'OFF')));
design = zeros(size(sensorData,3),length(onConds));
for k = 1:length(onConds)
    design(epoch_conditions(:,2)==onConds(k),k)=1;
end
