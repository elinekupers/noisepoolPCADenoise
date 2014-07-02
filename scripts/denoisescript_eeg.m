clear all;

%% load data
% get data into [channel x time x epoch] format 
% create corresponding design matrix [epoch x n] format, here n = 1
sessionnum  = 16;
conditionNumbers = 1:4;%3:4;
T = 1;
% Path to session, conditions
[sessionDir, conditionNames, conditionNumbers] = eegGetDataPaths(sessionnum, conditionNumbers);
% Paths to trials, headers
[pth, files, Axx, subj, thedate] = eegGetTrials(sessionDir, conditionNumbers);
% Extract trial info from a sample file
hdr = eegGetTrialHeader(pth, files{1}(1), T);
% Extract time series
[sensorData, dataTrials] = eegLoadtrials(pth, files, 'all', 'all', false, T);
% Remove bad epochs (epochs with more than 50% bad channels)
[sensorData,okEpochs] = eegRemoveBadEpochs(sensorData, 0.5);
% Replace bad channels (channels with more than 50% bad epochs)
badChannels = eegIdenitfyBadChannels(sensorData, 0.5);
fprintf('%d badChannels\n', sum(badChannels));
sensorData = eegReplaceBadChannels(sensorData, [], badChannels);
% Replace remaining bad epochs
sensorData = eegReplaceBadEpochs(sensorData);
% Check that we have eliminated all NaNs
disp(cellfun(@(x) mean(isnan(x(:))), sensorData))
% concatenate across conditions
tepochs    = [];
for ii = 1:length(conditionNames)
    tepochs = cat(1,tepochs,ii*ones(size(sensorData{ii},3),1));
end
sensorData = catcell(3,sensorData);
% design matrix
onConds = find(cellfun(@isempty,strfind(conditionNames,'off')));
design  = zeros(size(sensorData,3),length(onConds));
for k = 1:length(onConds)
    design(tepochs==onConds(k),k) = 1;
end
% save file directory 
eegDataDir = fileparts(pth);

% save(fullfile('/arc/1.2/p3/helena/denoisesuite/eegfigs/tmpeeg/',...
%     sprintf('%02d_%s',sessionnum,sessionDir)),'sensorData', 'design');
%clear allResults allEval
%size(sensorData)
%% Denoise 
% define some parameters for doing denoising 
fmax = 150; f = hdr.f;
%freq = eegGetSLandABfrequencies(f(f<fmax), max(hdr.t), 18);
freq = eegGetSLandABfrequencies((0:fmax), 1, 18);
evokedfun = @(x)getstimlocked(x,freq);
evalfun   = {@(x)getbroadband(x,freq), @(x)getstimlocked(x,freq)};
%evalfun   = {@(x)getbroadband(x,freq)};

npool = 60;
clear opt
opt.freq = freq;
opt.npcs = 45;
opt.xvalratio = -1;
opt.resampling = {'xval','xval'};
opt.npoolmethod = {'r2','n',npool};
%opt.npoolmethod = {'r2','thres',0};
opt.pccontrolmode = 0;

%opt.pcstop = -10;
opt.verbose = true;
% do denoising 
% use evokedfun to do noise pool selection 
% use evalfun   to do evaluation 
for ii = 0%:4 % loop through the different types of nulls 
    opt.pccontrolmode = ii;
    tic
    [results,evalout] = denoisedata(design,sensorData,evokedfun,evalfun,opt);
    %allResults{ii+1} = results;
    %allEval{ii+1} = evalout;
    toc    
end

%save(fullfile('tmpeeg',sprintf('%02d_%s_n%d',sessionnum,sessionDir,npool)),'results');

% doSave = true;
% if doSave
%     savename = fullfile('tmpeeg',sprintf('%s_n%d_nulls',sessionDir,npool));
%     fprintf('data saved: %s\n', savename);
%     save(savename,'allResults','allEval');
% end

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Do some evaluations 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% look at whether broadband signal as a function of pcs
warning off
printFigsToFile = true;
types = {'BroadBand','StimulusLocked'};
noisepool = results.noisepool;
opt = results.opt;

if printFigsToFile
    savepth = fullfile(eegDataDir, 'denoisefigures0');
    fprintf('Saving images to %s\n', savepth);
    if ~exist(savepth, 'dir'), mkdir(savepth); end
    if length(conditionNames)==2, stradd0 = conditionNames{1}; else stradd0 = 'ALL'; end
    stradd0 = sprintf('%s_n%d',stradd0,npool);
    %if opt.pccontrolmode, stradd0 = sprintf('%s_NULL%d',stradd0,opt.pccontrolmode); end
    disp(stradd0)
end

%%
% compute axses 
clims    = getblims(evalout,[15,85;5,95]);
numconds = size(clims,1);

if ~printFigsToFile
    whichbeta = 1;
    clims_ab = squeeze(clims(whichbeta,:,:))';
    
    figure('Position',[1 200 1200 800]);
    for p = 0:opt.npcs
        for fh = 1:length(evalfun)
            % [1 x channel x perms], averaged across perms  
            beta = mean(evalout(p+1,fh).beta(whichbeta,:,:),3); 
            r2   = evalout(p+1,fh).r2;
            
            subplot(2,2,fh);   % plot beta weights 
            ttl = sprintf('%s: PC = %02d', types{fh}, p);
            eegPlotMap(beta,[],'jet',ttl,'zbuffer',clims_ab(fh,:));
            
            subplot(2,2,fh+2); % plot r^2 
            ttl = sprintf('%s R2: PC = %02d', types{fh}, p);
            eegPlotMap(r2,[],'jet',ttl,'zbuffer',[]);
        end
        pause;
    end
    
else
    for p = 0:opt.npcs
        for fh = 1:length(evalfun)
            % loop through each beta and plot beta weights
            for whichbeta = 1:numconds
                
                beta = mean(evalout(p+1,fh).beta(whichbeta,:,:),3); % [1 x channel x perms], averaged across perms
                
                eegPlotMap(beta, 1,'jet',sprintf('%s: PC = %02d', types{fh}, p),'zbuffer', clims(whichbeta,:,fh));
                
                if numconds > 1, stradd = [stradd0, '_', conditionNames{onConds(whichbeta)}];
                else stradd = stradd0; end

                %saveas(fH,fullfile(savepth, sprintf('%s_%s_PC%02d.png',stradd,types{fh},p)),'png');
                figurewrite(sprintf('%s_%s_PC%02d',stradd,types{fh},p),[],[],savepth,0);
                
            end
            
            if fh == 1
                % plot r^2
                r2   = evalout(p+1,fh).r2;
                eegPlotMap(r2,1,'jet',sprintf('%s R2: PC = %02d', types{fh}, p),'zbuffer',[]);
                
                figurewrite(sprintf('%s_R2_%s_PC%02d',stradd0,types{fh},p),[],[],savepth,0);
            end
        end
    end
end

%% look at coverage of noise channels

signalnoiseVec = zeros(1,128);
signalnoiseVec(noisepool)  = 1;
eegPlotMap(signalnoiseVec,1,'autumn',sprintf('Noise channels: N = %d',sum(noisepool)),'zbuffer',[0,1]);
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
    vline(results.pcnum(fh),'k');
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
% this assumes we've computed this above 
%load(savename);
r2 = [];
for nn = 1:5
    evalout = allEval{nn};
    r2 = cat(3, r2, cat(1,evalout(:,1).r2)); 
end
noisepool = allResults{1}.noisepool;

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
    if kk == length(funcs), legend(nulltypes,'location','best'); end
end

if printFigsToFile
    figurewrite(sprintf('%s_Comparisons_R2vPCs_%s',stradd0,types{1}),[],[],savepth);
end

%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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