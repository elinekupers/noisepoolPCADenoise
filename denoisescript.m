
%% load data
% get data into [channel x time x epoch] format 
% create corresponding design matrix [epoch x n] format, here n = 1
megDataDir = '/Users/Helena/Work/EEG/MEG/data/05_SSMEG_04_04_2014/';
conditionNames = {'ON FULL','OFF FULL'};
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
evalfun   = {@(x)getbroadband(x,freq), @(x)getstimlocked(x,freq)};

opt.freq = freq;
opt.npcs = 50;
opt.xvalratio = -1;
opt.xvalmaxperm = 100;
opt.resampling = {'','xval'};
opt.npoolmethod = {'r2',[],'thres',0};
% do denoising 
% use evokedfun to do noise pool selection 
% use evalfun   to do evaluation 
[evalout,noisepool] = denoisedata(design,sensorData,evokedfun,evalfun,opt);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Do some evaluations 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% look at whether broadband signal as a function of pcs
printFigsToFile = true;
savepth = fullfile(megDataDir, 'denoisefigures0');
if printFigsToFile
    fprintf('Saving images to %s\n', savepth);
    if ~exist(savepth, 'dir'), mkdir(savepth); end
    stradd = sprintf('meg_%s', conditionNames{1}(:,4:end));
end

clims_ab = zeros(length(evalfun),2);
types = {'BroadBand','StimulusLocked'};
if ~printFigsToFile, figure('Position',[1 200 1200 600]); end
for p = 0:opt.npcs
    for fh = 1:2
        this_val = to157chan(mean(evalout(p+1,fh).beta),~badChannels, 'nans');
        
        if ~printFigsToFile, subplot(1,2,fh); cla; else clf; end
        ttl = sprintf('%s: PC = %02d', types{fh}, p);
        fH = megPlotMap(this_val,[],[],'jet',ttl);
        %ssm_plotOnMesh(this_val, [], [], [],'2d');
        
        % set axes 
        cl = get(gca, 'CLim');
        if p == 0, clims_ab(fh,:) = [-1, 1]*max(abs(cl))*1.2; end
        set(gca, 'CLim', clims_ab(fh,:));
        
        if printFigsToFile
            %saveas(fH,fullfile(savepth, sprintf('%s_%s_PC%02d.png',stradd,types{fh},p)),'png');
            figurewrite(sprintf('%s_%s_PC%02d.png',stradd,types{fh},p),[],[],savepth,0);
        end
    end
    if ~printFigsToFile, pause; end
end

%% look at coverage of noise channels

signalnoiseVec = zeros(1,size(sensorData,1));
signalnoiseVec(noisepool)  = 1;
signalnoiseVec = to157chan(signalnoiseVec,~badChannels,'nans');

fH = megPlotMap(signalnoiseVec,[-0.2,1],[],'autumn','Noise channels');
colorbar off;
if printFigsToFile, figurewrite('noisepool',[],[],savepth); end

%%
% look at how the r^2 changes as a function of denoising 
r2 = []; % npcs x channels [x evalfuns]
for fh = 1:length(evalfun)
    r2 = cat(3, r2,cat(1,evalout(:,fh).r2));
end

figure('Position',[1 200 600 600]);
fh = 1;
ax(1) = subplot(2,2,[1,2]);
imagesc(r2(:,:,fh)'); colorbar;
xlabel('n pcs'); ylabel('channel number');
title('R^2 as a function of denoising');

ax(2) = subplot(2,2,3);
plot(0:opt.npcs, r2(:,:,fh),'color',[0.5,0.5,0.5]); hold on;
plot(0:opt.npcs, mean(r2(:,:,fh),2),'r'); hold on;
xlabel('n pcs'); ylabel('r2');
title('R^2 for individual channels')

ax(3) = subplot(2,2,4);
plot(0:opt.npcs, mean(r2(:,:,fh),2),'b'); hold on;
plot(0:opt.npcs, mean(r2(:,~noisepool,fh),2),'r');
xlabel('n pcs'); ylabel('average r2');
legend('all channels','non-noise channels','Location','best');
title('mean R^2')

if printFigsToFile
    fs = 14;
    for ii = 1:3
        set(get(ax(ii),'Title'),'FontSize',fs);
        set(get(ax(ii),'XLabel'),'FontSize',fs);
        set(get(ax(ii),'YLabel'),'FontSize',fs);
        set(ax(ii),'box','off','tickdir','out','ticklength',[0.025 0.025]);
    end
    figurewrite('R2vPCs',[],[],savepth);
end

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