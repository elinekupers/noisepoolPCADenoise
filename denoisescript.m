
%% load data
% get data into [channel x time x epoch] format 
% create corresponding design matrix [epoch x n] format, here n = 1
megDataDir = '/Users/Helena/Work/EEG/MEG/data/05_SSMEG_04_04_2014/';
conditionNames = {'ON LEFT','OFF LEFT'};
numepochs  = zeros(1,2);
sensorData = [];
for ii = 1:2
    %dataName = ['ts_', lower(regexprep(conditionNames{ii},' ','_')), '_epoched'];
    dataName = ['ts_', lower(regexprep(conditionNames{ii},' ','_'))];
    disp(dataName);
    data     = load(fullfile(megDataDir,dataName));
    currdata = data.(dataName);
    numepochs(ii) = size(currdata,2);
    sensorData = cat(2,sensorData,currdata);
end
% format data into the right dimensions 
sensorData = permute(sensorData,[3,1,2]);
sensorData = sensorData(1:157,:,:);
% remove bad epochs
[sensorData,okEpochs] = megRemoveBadEpochs({sensorData},0.5);
% find bad channels 
badChannels = megIdenitfyBadChannels(sensorData, 0.5);
badChannels(98) = 1; % for now we add this in manually
fprintf('badChannels : %g \n', find(badChannels)');
% remove bad epochs and channels
net=load('meg160xyz.mat');
sensorData = megReplaceBadEpochs(sensorData,net);
sensorData = sensorData{1};
sensorData = sensorData(~badChannels,:,:);
% design matrix
design = zeros(size(sensorData,3),1);
design(1:numepochs(1),1) = 1;

%% Denoise 
% define some parameters for doing denoising 
T = 1; fmax = 150;
freq = megGetSLandABfrequencies((0:fmax)/T, T, 12/T);
evokedfun = @(x)getstimlocked(x,freq);
evalfun   = @(x)getbroadband(x,freq);

opt.freq = freq;
opt.npcs = 30;
opt.xvalratio = 0.99;
opt.xvalmaxperm = 100;
opt.resampling = {'','xval'};
% do denoising 
% use evokedfun to do noise pool selection 
% use evalfun   to do evaluation 
evalout2 = denoisedata(design,sensorData,evokedfun,evalfun,opt);

%%
% look at whether broadband signal at a function of pcs
for p = 0:opt.npcs
    this_ab = to157chan(mean(evalout(p+1).beta,1),~badChannels, 'nans');
    fH(3) = figure(3); clf;
    ssm_plotOnMesh(this_ab, [], [], [],'2d');
    colorbar; set(fH(3), 'renderer', 'zbuffer'); colormap jet;
    tH = title(sprintf('Broadband (PC = %02d)', p));
    cl = get(gca, 'CLim');
    if p == 0, clims_ab = [-1, 1]*max(abs(cl)); end
    set(gca, 'CLim', clims_ab);
    pause;
end

%%
% look at how the r^2 changes as a function of denoising 
r2 = cat(1,evalout.r2);
figure;
subplot(2,1,1);
imagesc(r2'); colorbar;
xlabel('n pcs'); ylabel('channel number');
subplot(2,1,2);
plot(0:opt.npcs, r2,'color',[0.5,0.5,0.5]); hold on;
plot(0:opt.npcs, mean(r2,2),'r'); hold on;
xlabel('n pcs'); ylabel('r2');

%%
r2o = cat(1,evalout0.r2);
r2  = cat(1,evalout2.r2);
axismin = -40; axismax = 50;
for p = 0:opt.npcs
    
    subplot(1,3,1)
    plot(r2o(p+1,:),r2(p+1,:),'or');
    line([axismin,axismax],[axismin,axismax],'color','k');
    xlim([axismin,axismax]); ylim([axismin,axismax]); axis square;
    
    subplot(1,3,[2,3]); hold off; 
    %plot(r2o(p+1,:)-r2(p+1,:));
    plot(r2o(p+1,:),'b'); hold on;
    plot(r2(p+1,:),'r');
    ylim([axismin,axismax]);
    
    title(sprintf('PC = %d',p));
    pause;
end