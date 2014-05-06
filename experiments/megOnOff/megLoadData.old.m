function [sensorData, badChannels, tepochs, epochGroup] = ...
    megLoadData(megDataDir,conditionNames,opt)
% load meg data from data and format appropriately 

if notDefined('opt'),    opt = struct(); end
if ~isfield(opt,'alt_epoch'),    opt.alt_epoch = 0;  end
if ~isfield(opt,'group_epoch'),  opt.group_epoch  = 1;  end
if ~isfield(opt,'shuffle_epoch'), opt.shuffle_epoch = false; end
if ~isfield(opt,'remove_strtend_epoch'), opt.remove_strtend_epoch = false; end
if ~isfield(opt,'badepoch_avgchannum'),  opt.badepoch_avgchannum = 6; end
if ~isfield(opt,'verbose'),      opt.verbose = true; end

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

% format data into the right dimensions
sensorData = permute(sensorData,[3,1,2]);
sensorData = sensorData(1:157,:,:);
% replace missing data with nan's
sensorData(sensorData==0) = nan; 

% switch into alternating conditions %<-- this needs to be updated
if opt.alt_epoch > 0
    if opt.verbose, fprintf('(megLoadData) alternating every %d epochs\n', opt.alt_epoch); end
    epochinds = reshape([1:360]',opt.alt_epoch,[],2);   % group 6 epochs at a time
    epochinds = permute(epochinds,[1,3,2]); % 6 ON, 6 OFF
    epochinds = epochinds(:);
    sensorData = sensorData(:,:,epochinds);
    tepochs = tepochs(epochinds);
end

% group by trial
if opt.group_epoch > 1
    if opt.verbose, fprintf('(megLoadData) grouping everying %d epochs for denoising\n', opt.group_epoch); end
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
    sensorData = sensorData(:,:,kept_epochs(:));
    tepochs = tepochs(kept_epochs(:));
    epochGroup = epochGroup(kept_epochs(:));
end

% remove bad epochs
[sensorData,okEpochs] = megRemoveBadEpochs({sensorData},0.5);
tepochs = tepochs(okEpochs{1});
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
