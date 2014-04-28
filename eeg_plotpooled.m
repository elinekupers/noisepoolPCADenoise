clear all;
sessionNums = 8:18;
rootDir = strrep(which('setup.m'),'denoisesuite/setup.m','');
rootDir = fullfile(rootDir,'EEG','Data','PowerDiva');

%% check for files 
filestem = '*Comparisons_R2vPCs_BroadBand*';
for k = 1:length(sessionNums)
    fprintf(' session %d \n', sessionNums(k));
    [sessionDir, conditionNames, conditionNumbers] = eegGetDataPaths(sessionNums(k));
    denoiseDir = fullfile(rootDir,sessionDir,'denoisefigures0');
    dir(fullfile(denoiseDir,filestem))
    fprintf('====================\n\n');
end

%% plot R^2 as a function of pcs

%figure('position',[1,600,800,300]); 
figure('position',[1,600,800,600]); 
fs = 16;
plotType = 3;
printFigsToFile = false;  
optpcs = zeros(1,length(sessionNums));
allPCchan = [];

for k = 1:length(sessionNums)
    fprintf(' session %d \n', sessionNums(k));
    [sessionDir, conditionNames, conditionNumbers] = eegGetDataPaths(sessionNums(k));
    filename = fullfile('tmpeeg',sprintf('%s*_nulls*',sessionDir));
    files = dir(filename);
    
    thisfile = fullfile('tmpeeg',files(end).name); disp(thisfile);
    load(thisfile);
    fprintf(' done loading\n');
    
    % look at r2 as a function of pcs 
    %-------------------------------------------
    xvaltrend = [];
    for pcnull = 0:4
        results = allResults{pcnull+1};
        evalout = allEval{pcnull+1};
        opt = results.opt; noisepool = results.noisepool;
        if pcnull ==0, disp(opt.npoolmethod), end
        
        r2 = []; % npcs x channels [x evalfuns]
        for fh = 1%:length(evalfun)
            r2 = cat(3, r2,cat(1,evalout(:,fh).r2));
        end
        if pcnull == 0
            % top x number of pcs
            pcchan = false(size(noisepool));
            maxr2 = max(r2,[],1); % max cross validation for each channel
            [~, idx] = sort(maxr2,'descend');
            pcchan(idx(1:min(10,length(idx)))) = 1;
        end
        xvaltrend = cat(2, xvaltrend, mean(r2(:,pcchan,1),2));
    end
    
    AllPCtrend{k} = xvaltrend;  % aggregates here 
    allPCchan = cat(1,allPCchan,pcchan);
    
    chosen = choosepc(xvaltrend(:,1),1.05);
    optpcs(k) = chosen;
    
    % plot
    if plotType == 1
        ax(1) = subplot(1,3,1); cla;
        plot(0:opt.npcs, r2(:,:,1));
        title('R^2 for individual channels')
        
        ax(2) = subplot(1,3,2); cla;
        plot(0:opt.npcs, mean(r2(:,:,1),2),'b','linewidth',2); hold on;
        plot(0:opt.npcs, mean(r2(:,~noisepool,1),2),'r','linewidth',2);
        title('mean R^2: all; non-noise');
        
        ax(3) = subplot(1,3,3); cla
        plot(0:opt.npcs, xvaltrend(:,1), 'k','linewidth',2);
        title('mean R^2: top 10')
        
        for ii = 1:3
            subplot(1,3,ii);
            xlabel('n pcs'); ylabel('R2');
            vline(results.pcnum(1),'r');
            vline(chosen,'k');
            xlim([0,50]);
            axis square;
            makeprettyaxes(gca,12);
        end
        suptitle(sprintf('EEG Session %d : %s', sessionNums(k), sessionDir));
        
    elseif plotType == 3
        
        subplot(3,4,k); cla;
        plot(0:opt.npcs,r2(:,:,1));
        title(sprintf('s%d : %s', sessionNums(k), sessionDir));
        vline(chosen,'k');
        xlim([0,50]); axis square;
        makeprettyaxes(gca,12);
    end
    
    %pause;
    clear allEval
    fprintf('====================\n\n');
    
    if printFigsToFile
        if plotType == 1
            figurewrite(sprintf('%02d_%s', sessionNums(k), sessionDir),[],[],'eegfigs',1);
        end
    end    
end
% save optpcs optpcs sessionNums
if printFigsToFile && plotType == 3
    figurewrite('PCselection_allsubjs2',[],[],'eegfigs',1);
end

%%
% plot R^2 as a function of PCs, with all subjects together
% also plot the different types of nulls 

trends = catcell(3,AllPCtrend);
fudge = [0,3,3,3,3];
colors = copper(12);
ttls = {'Original','Phase scrambled','Order shuffled','Amplitude scrambled','Random pcs'};
for k = 1:5
    subplot(2,4,k+fudge(k)); cla; hold on;
    for nn = 1:size(trends,3)-1 % loop through subjects
        plot(0:45,squeeze(trends(:,k,nn)),'color',colors(nn,:));
    end
    xlabel('n pcs'); ylabel('R2'); xlim([0,50]); ylim([-5,25]); axis square;
    title(ttls{k}); makeprettyaxes;
    
end
if printFigsToFile
    figurewrite('PCselection_allsubjs',[],[],'eegfigs',1);
end

%% plot the spatial maps as a function of denoising 
%figure('position',[1,600,1600 500]); fs = 16;
printFigsToFile = true;  
optpcs = zeros(1,length(sessionNums));
plotSNR = true;

allSNR1 = []; allSNR2 = []; allNoisepool = [];
plotType = 'noise';
for k = 1:length(sessionNums)
    fprintf(' session %d \n', sessionNums(k));
    [sessionDir, conditionNames, conditionNumbers] = eegGetDataPaths(sessionNums(k));
    filename = fullfile('tmpeeg',sprintf('%02d_%s*',sessionNums(k),sessionDir)); 
    files = dir(filename);
    
    thisfile = fullfile('tmpeeg',files(end).name); disp(thisfile);
    load(thisfile); fprintf(' done loading\n');

    opt = results.opt; noisepool = results.noisepool;
    disp(opt.npoolmethod);
    %disp(results.pcnum(1))
    
    if plotSNR
        sl_snr1 = abs(results.origmodel(2).beta_md)./results.origmodel(2).beta_se;
        ab_snr1 = abs(results.origmodel(1).beta_md)./results.origmodel(1).beta_se;
        ab_snr2 = abs(results.finalmodel(1).beta_md)./results.finalmodel(1).beta_se;
        clims_sl = [0, max(sl_snr1)];
        clims_ab = [0, max([ab_snr1, ab_snr2])];
        
    else
        sl_snr1 = results.origmodel(2).r2;
        ab_snr1 = results.origmodel(1).r2;
        ab_snr2 = results.finalmodel(1).r2;
        clims_sl = [min(sl_snr1), max(sl_snr1)];
        clims_ab = [min([ab_snr1, ab_snr2]), max([ab_snr1, ab_snr2])];
    end
    
    allSNR1 = cat(1,allSNR1,ab_snr1);
    allSNR2 = cat(1,allSNR2,ab_snr2);
    allNoisepool = cat(1,allNoisepool, noisepool);
    
    switch plotType
        case 'SNR'
            subplot(1,3,1);
            eegPlotMap(sl_snr1,[],'jet','Stimulus Locked Original','zbuffer',clims_sl);
            subplot(1,3,2);
            eegPlotMap(ab_snr1,[],'jet','Broad Band Original','zbuffer',clims_ab);
            subplot(1,3,3);
            eegPlotMap(ab_snr2,[],'jet',sprintf('Broad Band PC %d',results.pcnum(1)),'zbuffer',clims_ab);
            
            %suptitle(sprintf('EEG Session %d : %s', sessionNums(k), sessionDir));
        case 'noise'
            fH = eegPlotMap(double(noisepool),[],'autumn',sprintf('%s: N=%d',sessionDir, sum(noisepool)),[],[0,1]);
            colorbar off
    end
            
    %pause;
    %clear results;
    fprintf('====================\n\n');
    
    if printFigsToFile
        if strcmp(plotType,'SNR')
            if plotSNR, header = 'SNR'; else header = 'R2'; end
            figurewrite(sprintf('%sMap_%02d_%s', header, sessionNums(k), sessionDir),[],[],'eegfigs',1);
        elseif strcmp(plotType,'noise')
            figurewrite(sprintf('noisepool%d_%s',sessionNums(k),sessionDir),[],[],'eegfigs',1);
        end
    end 
end

%% Plot SNR before and after 

axismin = 0; axismax = 20;
colors = jet(12);
allNoisepool = logical(allNoisepool);
hold on;
for nn = 1:length(sessionNums)-1
    fprintf(' session %d \n', sessionNums(nn));
    sessionDir = eegGetDataPaths(sessionNums(nn));
    
    subplot(3,4,nn); cla;
    plot(allSNR1(nn,:),allSNR2(nn,:),'o','color',[1,1,1]*0.2); hold on;
    plot(allSNR1(nn,allNoisepool(nn,:)),allSNR2(nn,allNoisepool(nn,:)),'or');
    axismax = max([allSNR1(nn,:),allSNR2(nn,:)])*1.2;
    line([axismin,axismax],[axismin,axismax],'color','k');
    xlim([axismin,axismax]); ylim([axismin,axismax]); axis square;
    xlabel('orig model SNR'); ylabel('final model SNR');
    title(sprintf('S%d : %s', sessionNums(nn), sessionDir))
    makeprettyaxes(gca,12);
end

if printFigsToFile
    figurewrite('SNRbeforeafter_allsubjs',[],[],'eegfigs',1);
end

%% SNR before and after, mean across good PCs
% depends on having computed allPCchan
allPCchan = logical(allPCchan);
for nn = 1:length(sessionNums)-1
    mean_before(nn) = mean(allSNR1(nn,allPCchan(nn,:)),2);
    mean_after(nn) = mean(allSNR2(nn,allPCchan(nn,:)),2);
end
colors = copper(12);
tmp = [mean_before; mean_after]';
hold on;
for nn = 1:length(sessionNums)-1
    plot(1:2,tmp(nn,:),'o-','color',colors(nn,:));
end
xlim([0,3]);
set(gca,'xtick',1:2,'xticklabel',{'Before','After'}); 
ylabel('SNR');
makeprettyaxes(gca,14); axis square;

if printFigsToFile
    figurewrite('SNRbeforeafter_allsubjs2',[],[],'eegfigs',1);
end