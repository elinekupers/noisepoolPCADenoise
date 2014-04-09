
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
evalfun   = {@(x)getbroadband(x,freq), @(x)getbroadbandlog(x,freq)};

opt.freq = freq;
opt.npcs = 50;
opt.xvalratio = -1;
opt.xvalmaxperm = 100;
opt.resampling = {'','xval'};
% do denoising 
% use evokedfun to do noise pool selection 
% use evalfun   to do evaluation 
[evalout,noisepool] = denoisedata(design,sensorData,evokedfun,evalfun,opt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Do some evaluations 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% look at whether broadband signal at a function of pcs
figure('Position',[1 200 1200 600]);
clims_ab = zeros(length(evalfun),2);
for p = 0:opt.npcs
    for fh = 1:2
        subplot(1,2,fh); cla; 
        this_ab = to157chan(mean(evalout(p+1,fh).beta),~badChannels, 'nans');
        
        ssm_plotOnMesh(this_ab, [], [], [],'2d');
        
        colorbar% off; 
        colormap jet;
        tH = title(sprintf('Broadband (PC = %02d)', p));
        cl = get(gca, 'CLim');
        if p == 0, clims_ab(fh,:) = [-1, 1]*max(abs(cl)); end
        set(gca, 'CLim', clims_ab(fh,:));
    end
    pause;
end

%%
% look at how the r^2 changes as a function of denoising 
r2 = []; % npcs x channels [x evalfuns]
for fh = 1:length(evalfun)
    r2 = cat(3, r2,cat(1,evalout(:,fh).r2));
end

figure('Position',[1 200 600 600]);
fh = 1;
subplot(2,2,[1,2]);
imagesc(r2(:,:,fh)'); colorbar;
xlabel('n pcs'); ylabel('channel number');
title('R^2 as a function of denoising');

subplot(2,2,3);
plot(0:opt.npcs, r2(:,:,fh),'color',[0.5,0.5,0.5]); hold on;
plot(0:opt.npcs, mean(r2(:,:,fh),2),'r'); hold on;
xlabel('n pcs'); ylabel('r2');

subplot(2,2,4);
plot(0:opt.npcs, mean(r2(:,:,fh),2),'b'); hold on;
plot(0:opt.npcs, mean(r2(:,~noisepool,fh),2),'r');
xlabel('n pcs'); ylabel('average r2');
legend('all channels','non-noise channels','Location','best');

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