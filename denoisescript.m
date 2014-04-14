clear all;
%% load data
% get data into [channel x time x epoch] format 
% create corresponding design matrix [epoch x n] format, here n = 1
rootDir = strrep(which('setup.m'),'denoisesuite/setup.m','');
dataset = '03_SSMEG_03_31_2014';
megDataDir = fullfile(rootDir,'data',dataset);
conditionNames = {'ON FULL','OFF FULL','ON LEFT','OFF LEFT','ON RIGHT','OFF RIGHT'};
%conditionNames = {'ON RIGHT','OFF RIGHT'};
tepochs    = [];
sensorData = [];
for ii = 1:length(conditionNames)
    %dataName = ['ts_', lower(regexprep(conditionNames{ii},' ','_')), '_epoched'];
    dataName = ['ts_', lower(regexprep(conditionNames{ii},' ','_'))];
    disp(dataName);
    data     = load(fullfile(megDataDir,dataName));
    currdata = data.(dataName);
    tepochs  = cat(1,tepochs,ii*ones(size(currdata,2),1));
    sensorData = cat(2,sensorData,currdata);
end
% format data into the right dimensions 
sensorData = permute(sensorData,[3,1,2]);
sensorData = sensorData(1:157,:,:);
% remove bad epochs
[sensorData,okEpochs] = megRemoveBadEpochs({sensorData},0.5);
tepochs = tepochs(okEpochs{1});
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
onConds = find(cellfun(@isempty,strfind(conditionNames,'OFF')));
design = zeros(size(sensorData,3),length(onConds));
for k = 1:length(onConds)
    design(tepochs==onConds(k),k) = 1;
end

%save(sprintf('tmp/%s',dataset),'sensorData', 'design', 'freq', 'badChannels');

%% Denoise 
% define some parameters for doing denoising 
T = 1; fmax = 150;
freq = megGetSLandABfrequencies((0:fmax)/T, T, 12/T);
evokedfun = @(x)getstimlocked(x,freq);
evalfun   = {@(x)getbroadband(x,freq), @(x)getstimlocked(x,freq)};
%evalfun   = @(x)getbroadband(x,freq);

opt.freq = freq;
opt.npcs = 50;
opt.xvalratio = -1;
opt.xvalmaxperm = 100;
opt.resampling = {'','xval'};
%opt.npoolmethod = {'r2',[],'n',60};
opt.npoolmethod = {'r2',[],'thres',0};
% do denoising 
% use evokedfun to do noise pool selection 
% use evalfun   to do evaluation 
[finalmodel,evalout,noisepool,denoisedspec] = denoisedata(design,sensorData,evokedfun,evalfun,opt);

%return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Do some evaluations 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% look at whether broadband signal as a function of pcs
warning off
printFigsToFile = true;

if printFigsToFile
    savepth = fullfile(megDataDir, 'denoisefigures0');
    fprintf('Saving images to %s\n', savepth);
    if ~exist(savepth, 'dir'), mkdir(savepth); end
end

% compute axses 
clims = getblims(evalout,[15,85;5,95]);

types = {'BroadBand','StimulusLocked'};
if ~printFigsToFile
    whichbeta = 1;
    clims_ab = squeeze(clims(whichbeta,:,:))';
    
    figure('Position',[1 200 1200 800]);
    for p = 0:opt.npcs
        for fh = 1:length(evalfun)
            beta = mean(evalout(p+1,fh).beta(whichbeta,:,:),3); % [1 x channel x perms], averaged across perms  
            beta = to157chan(beta,~badChannels, 'nans');          % map back to 157 channel space
            r2   = evalout(p+1,fh).r2;
            r2   = to157chan(r2,~badChannels,'nans');
            
            subplot(2,2,fh);   % plot beta weights 
            ttl = sprintf('%s: PC = %02d', types{fh}, p);
            fH = megPlotMap(beta,clims_ab(fh,:),[],'jet',ttl);
            %ssm_plotOnMesh(beta, [], [], [],'2d');
            
            subplot(2,2,fh+2); % plot r^2 
            ttl = sprintf('%s R2: PC = %02d', types{fh}, p);
            fH = megPlotMap(r2,[],[],'jet',ttl);
        end
        pause;
    end
    
else
    for p = 0:opt.npcs
        for fh = 1:length(evalfun)
            % loop through each beta and plot beta weights
%             for whichbeta = 1:3
%                 
%                 beta = mean(evalout(p+1,fh).beta(whichbeta,:,:),3); % [1 x channel x perms], averaged across perms
%                 beta = to157chan(beta,~badChannels, 'nans');        % map back to 157 channel space
%                 
%                 ttl = sprintf('%s: PC = %02d', types{fh}, p);
%                 fH = megPlotMap(beta,clims(whichbeta,:,fh),[],'jet',ttl);
%                 
%                 stradd = sprintf('ALL_%s', conditionNames{onConds(whichbeta)}(:,4:end));
%                 %saveas(fH,fullfile(savepth, sprintf('%s_%s_PC%02d.png',stradd,types{fh},p)),'png');
%                 figurewrite(sprintf('%s_%s_PC%02d',stradd,types{fh},p),[],[],savepth,0);
%                 
%             end
            
            if fh == 1
                % plot r^2
                r2   = evalout(p+1,fh).r2;
                r2   = to157chan(r2,~badChannels,'nans');
                ttl = sprintf('%s R2: PC = %02d', types{fh}, p);
                fH = megPlotMap(r2,[],[],'jet',ttl);
                
                figurewrite(sprintf('ALL_R2_%s_PC%02d',types{fh},p),[],[],savepth,0);
            end
        end
    end
    
end

%% look at coverage of noise channels

signalnoiseVec = zeros(1,size(beta,2));
signalnoiseVec(noisepool)  = 1;
signalnoiseVec = to157chan(signalnoiseVec,~badChannels,'nans');

fH = megPlotMap(signalnoiseVec,[0,1],[],'autumn',sprintf('Noise channels: N = %d',sum(noisepool)));
colorbar off;

if length(conditionNames)==2, stradd = conditionNames{1}(:,4:end); else stradd = 'ALL'; end
if printFigsToFile, figurewrite(sprintf('%s_noisepool',stradd),[],[],savepth); end


%% look at how the r^2 changes as a function of denoising 
r2 = []; % npcs x channels [x evalfuns]
for fh = 1:length(evalfun)
    r2 = cat(3, r2,cat(1,evalout(:,fh).r2));
end

for fh = 1:length(evalfun)
    figure('Position',[1 200 600 600]);
    
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
    vline(finalmodel(fh).pcnum,'k');
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
        figurewrite(sprintf('%s_R2vPCs_%s',stradd,types{fh}),[],[],savepth);
    else
        pause;
    end
end

return;
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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