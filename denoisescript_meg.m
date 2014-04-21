clear all;
%% load data
% get data into [channel x time x epoch] format 
% create corresponding design matrix [epoch x n] format, here n = 1
rootDir = strrep(which('setup.m'),'denoisesuite/setup.m','');
dataset = '02_SSMEG_02_28_2014';
megDataDir = fullfile(rootDir,'data',dataset);
conditionNames = {'ON FULL','OFF FULL','ON LEFT','OFF LEFT','ON RIGHT','OFF RIGHT'};
%conditionNames = {'ON FULL','OFF FULL'};
tepochs    = [];
sensorData = [];
for ii = 1:length(conditionNames)
    dataName = ['ts_', lower(regexprep(conditionNames{ii},' ','_')), '_epoched'];
    %dataName = ['ts_', lower(regexprep(conditionNames{ii},' ','_'))];
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
opt.pccontrolmode = 4;
opt.verbose = true;
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
    if length(conditionNames)==2, stradd0 = conditionNames{1}(:,4:end); else stradd0 = 'ALL'; end
    if opt.pccontrolmode, stradd0 = sprintf('%s_NULL%d',stradd0,opt.pccontrolmode); end
    disp(stradd0);
end

%%
% compute axses 
clims    = getblims(evalout,[15,85;5,95]);
numconds = size(clims,1);

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
            for whichbeta = 1:numconds
                
                beta = mean(evalout(p+1,fh).beta(whichbeta,:,:),3); % [1 x channel x perms], averaged across perms
                beta = to157chan(beta,~badChannels, 'nans');        % map back to 157 channel space
                
                ttl = sprintf('%s: PC = %02d', types{fh}, p);
                fH = megPlotMap(beta,clims(whichbeta,:,fh),[],'jet',ttl);
                
                if numconds > 1, stradd = [stradd0, '_', conditionNames{onConds(whichbeta)}(:,4:end)];
                else stradd = stradd0; end
                
                %saveas(fH,fullfile(savepth, sprintf('%s_%s_PC%02d.png',stradd,types{fh},p)),'png');
                figurewrite(sprintf('%s_%s_PC%02d',stradd,types{fh},p),[],[],savepth,0);
                
            end
            
            if fh == 1
                % plot r^2
                r2   = evalout(p+1,fh).r2;
                r2   = to157chan(r2,~badChannels,'nans');
                ttl = sprintf('%s R2: PC = %02d', types{fh}, p);
                fH = megPlotMap(r2,[],[],'jet',ttl);
                
                figurewrite(sprintf('%s_R2_%s_PC%02d',stradd0,types{fh},p),[],[],savepth,0);
            end
        end
    end
end

%% look at coverage of noise channels

signalnoiseVec = zeros(1,size(evalout(1,1).beta,2));
signalnoiseVec(noisepool)  = 1;
signalnoiseVec = to157chan(signalnoiseVec,~badChannels,'nans');

fH = megPlotMap(signalnoiseVec,[0,1],[],'autumn',sprintf('Noise channels: N = %d',sum(noisepool)));
colorbar off;

if printFigsToFile, figurewrite(sprintf('%s_noisepool',stradd0),[],[],savepth); end


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
        figurewrite(sprintf('%s_R2vPCs_%s',stradd0,types{fh}),[],[],savepth);
    else
        pause;
    end
end

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot the comparisons between the different kinds of null 

r2 = [];
for nn = 0:4
    if nn == 0
        filename = sprintf('tmp/%s_fitfull',dataset);
    else
        filename = sprintf('tmp/%s_fitfull_null%d',dataset,nn);
    end
    disp(filename); load(filename);
    r2 = cat(3, r2, cat(1,evalout(:,1).r2)); 
end

%%
figure('Position',[1 200 1000 500]);
npcs = size(r2,1)-1;
nulltypes = {'original','phase scrambled','order shuffled','amplitude scrambled','random pcs'};
for nn = 1:5
    subplot(2,4,nn); hold on;
    plot(0:npcs, r2(:,:,nn)); hold on;
    if nn == 1, ylims = get(gca,'ylim'); end
    ylim(ylims); xlim([0,50]); xlabel('npcs'); ylabel('R^2');
    axis square; title(nulltypes{nn});
end
colors = {'k','b','r','g','m'};
%top10 = r2(end,:,1) > prctile(r2(end,:,1),90,2);
ttls = {'mean(all)','mean(non-noise)'};
funcs = {@(x)mean(x,2), @(x)mean(x(:,~noisepool),2)}; %@(x)prctile(x(:,:,1),90,2)
for kk = 1:length(funcs)
    subplot(2,4,5+kk); hold on;
    for nn = 1:5
        curr_r = r2(:,:,nn);
        plot(0:npcs,funcs{kk}(curr_r),colors{nn},'linewidth',2);
    end
    axis square; title(ttls{kk});
    xlim([0,50]); xlabel('npcs'); ylabel('R^2');
    if kk == length(funcs), legend(nulltypes,'location','bestoutside'); end
end

if printFigsToFile
    figurewrite(sprintf('ALLComparisons_R2vPCs_%s',types{1}),[],[],savepth);
end


%% %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% r2o = cat(1,evalout0.r2);
% r2  = cat(1,evalout2.r2);
% axismin = -40; axismax = 50;
% for p = 0:opt.npcs
%     
%     subplot(1,3,1)
%     plot(r2o(p+1,:),r2(p+1,:),'or');
%     line([axismin,axismax],[axismin,axismax],'color','k');
%     xlim([axismin,axismax]); ylim([axismin,axismax]); axis square;
%     
%     subplot(1,3,[2,3]); hold off; 
%     %plot(r2o(p+1,:)-r2(p+1,:));
%     plot(r2o(p+1,:),'b'); hold on;
%     plot(r2(p+1,:),'r');
%     ylim([axismin,axismax]);
%     
%     title(sprintf('PC = %d',p));
%     pause;
% end