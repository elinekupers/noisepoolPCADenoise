clear all;
sessionNums = 1:4;
tmpmegdir = '/Volumes/HelenaBackup/denoisesuite/tmpmeg/';

%% plot R^2 as a function of pcs
% long loading time
fs = 16;
printFigsToFile = false;
optpcs = zeros(1,length(sessionNums));
plotType = 2;
loadNull = false;

if plotType == 1
    figure('position',[1,600,800,300]);
elseif plotType == 2
    figure('position',[1,600,400,800]);
end

for k = 1:length(sessionNums)
    fprintf(' session %d \n', sessionNums(k));
    [sessionDir,conditionNames,megDataDir] = megGetDataPaths(sessionNums(k), 1:6);
    thisfile = fullfile(tmpmegdir,sprintf('%s_fitfull',sessionDir));
    %thisfile = fullfile(tmpmegdir,sprintf('%s_fitfull75',sessionDir));
    disp(thisfile); load(thisfile);
    fprintf(' done loading\n');
    if exist('result','var')
        noisepool = results.noisepool;
        opt = results.opt;
    end
    %
    % look at r2 as a function of pcs
    %-------------------------------------------
    xvaltrend = [];
    r2 = cat(1,evalout(:,1).r2); % npcs x channels
    % top x number of pcs
    pcchan = false(size(noisepool));
    maxr2 = max(r2,[],1); % max cross validation for each channel
    [~, idx] = sort(maxr2,'descend');
    pcchan(idx(1:min(10,length(idx)))) = 1;
    xvaltrend = cat(2, xvaltrend, mean(r2(:,pcchan,1),2));
    chosen = choosepc(xvaltrend(:,1),1.05);
    
    if loadNull
        allR2 = r2;
        for ii = 1:4
            thisfile = fullfile(tmpmegdir,sprintf('%s_fitfull_null%d',sessionDir,ii)); 
            disp(thisfile); load(thisfile);
            allR2 = cat(3,allR2,cat(1,evalout(:,1).r2));
        end
        allR2null{k} = allR2;
    end
    % aggregates here
    AllPCtrend{k} = xvaltrend;
    allPCchan{k}  = pcchan;
    optpcs(k) = chosen; 
    disp(opt.npcs); opt.npcs = size(r2,1)-1;
    
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
            xlabel('n pcs'); ylabel('R2'); axis square;
            xlim([0,50]);
            %vline(results.pcnum(1),'r');
            vline(finalmodel(1).pcnum,'r');
            vline(chosen,'k');
            makeprettyaxes(gca,12);
        end
        suptitle(sprintf('MEG Session %d : %s', sessionNums(k), sessionDir));
        
    elseif plotType == 2
        
        subplot(2,1,1);
        plot(0:opt.npcs, r2(:,:,1));
        
        subplot(2,1,2);
        plot(0:opt.npcs, xvaltrend(:,1), 'k','linewidth',2);
        
        for ii = 1:2
            subplot(2,1,ii);
            %xlabel('n pcs'); ylabel('R2'); 
            axis square; 
            xlim([0,50]); %tmp = get(gca,'ylim'); ylim([-1,tmp(2)]);
            makeprettyaxes(gca,14,14);
        end
        vline(chosen,'k');
        suptitle(sprintf('S%d : %s', sessionNums(k), sessionDir));
    end
    
    %pause;
    clear allEval results noisepool opt 
    fprintf('====================\n\n');
    
    if printFigsToFile
        figurewrite(sprintf('%02dn75_%s', sessionNums(k), sessionDir),[],[],'megfigs',1);
    end
end

%%
% plot R^2 as a function of PCs, with all subjects together
% also plot the different types of nulls 
figure('position',[1,600,1000,400]); 
fudge = [0,3,3,3,3];
ttls = {'Original','Phase scrambled','Order shuffled','Amplitude scrambled','Random pcs'};
plotType = 2;

switch plotType
    case 1 % plot all subjs in each panel
        colors = copper(4);
        for nn = 1:4 % loop through subjects
            r2 = allR2null{nn};
            for k = 1:5 % conditions
                subplot(2,4,k+fudge(k)); hold on;
                plot(0:size(r2,1)-1, mean(r2(:,allPCchan{nn},k),2),'color',colors(nn,:));
                
                if nn == 4
                    xlabel('n pcs'); ylabel('R2'); xlim([0,50]); ylim([-5,25]); axis square;
                    title(ttls{k}); makeprettyaxes;
                end
            end
        end
        if printFigsToFile
            figurewrite('PCselection_allsubjs',[],[],'megfigs',1);
        end
    case 2 % plot each subject separately 
        colors = ['k','b','r','g','m'];
        for nn = 1:4 % for each subject 
            r2 = allR2null{nn}; % npcs x nchannels x 5 conditions 
            subplot(1,4,nn);
            cla; hold on;
            for k = 1:5 % conditions
                 plot(0:size(r2,1)-1, mean(r2(:,allPCchan{nn},k),2),colors(k));
            end
            
            axis square; xlim([0,50]); tmp = get(gca,'ylim'); ylim([-1,tmp(2)]);
            title(sprintf('S%d', sessionNums(nn)));
            makeprettyaxes(gca,14,12);
            if nn == 4, legend(ttls,'location','bestoutside'); end 
        end
        if printFigsToFile
            figurewrite('PCselection_allsubjs2',[],-1,'megfigs',1);
        end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot the spatial maps as a function of denoising
fs = 16;
printFigsToFile = false;
optpcs = zeros(1,length(sessionNums));
plotType = '';
plotSNR  = true;

if strcmp(plotType, 'SNR'), figure('position',[1,600,1600 500]);  end
allSNR1 = []; allSNR2 = []; allNoisepool = [];
whichbetas = 3;

condNames = {'FULL','LEFT','RIGHT'};
for k = 1:length(sessionNums)
    fprintf(' session %d \n', sessionNums(k));
    [sessionDir,conditionNames,megDataDir] = megGetDataPaths(sessionNums(k), 1:6);
    thisfile = fullfile('megfigs',sprintf('%02d_%s',sessionNums(k),sessionDir));
    disp(thisfile);
    load(thisfile); fprintf(' done loading\n');
    
    opt = results.opt; noisepool = results.noisepool;
    disp(opt.npoolmethod);
    
    if plotSNR
        sl_snr1 = max(abs(results.origmodel(2).beta_md(whichbetas,:)),[],1)...
            ./mean(results.origmodel(2).beta_se(whichbetas,:),1);
        ab_snr1 = max(abs(results.origmodel(1).beta_md(whichbetas,:)),[],1)./...
            mean(results.origmodel(1).beta_se(whichbetas,:),1);
        ab_snr2 = max(abs(results.finalmodel(1).beta_md(whichbetas,:)),[],1)./...
            mean(results.finalmodel(1).beta_se(whichbetas,:),1);
        clims_sl = [0, max(sl_snr1)];
        clims_ab = [0, max([ab_snr1, ab_snr2])];
        
    else
        sl_snr1 = results.origmodel(2).r2;
        ab_snr1 = results.origmodel(1).r2;
        ab_snr2 = results.finalmodel(1).r2;
        clims_sl = [min(sl_snr1), max(sl_snr1)];
        clims_ab = [min([ab_snr1, ab_snr2]), max([ab_snr1, ab_snr2])];
    end
    
    sl_snr1a = to157chan(sl_snr1,~badChannels,'nans');
    ab_snr1a = to157chan(ab_snr1,~badChannels,'nans');
    ab_snr2a = to157chan(ab_snr2,~badChannels,'nans');
    noise2   = to157chan(noisepool,~badChannels,0);
    allSNR1 = cat(1,allSNR1,ab_snr1a);
    allSNR2 = cat(1,allSNR2,ab_snr2a);
    allNoisepool = cat(1,allNoisepool, noise2);
    
    switch plotType
        case 'SNR'
            subplot(1,3,1);
            megPlotMap(sl_snr1a,clims_sl,[],'jet','Stimulus Locked Original');
            subplot(1,3,2);
            megPlotMap(ab_snr1a,clims_ab,[],'jet','Broad Band Original');
            subplot(1,3,3);
            megPlotMap(ab_snr2a,clims_ab,[],'jet',sprintf('Broad Band PC %d',results.pcnum(1)));
            
        case 'noise'
            fH = megPlotMap(noise2,[0,1],[],'autumn',sprintf('Noise channels: N = %d',sum(noisepool)));
            colorbar off;
    end
    
    %pause;
    clear results badChannels
    fprintf('====================\n\n');
    
    if printFigsToFile
        if strcmp(plotType,'SNR')
            if plotSNR, header = 'SNR'; else header = 'R2'; end
            fname = sprintf('%sMap_%02d_%s', header, sessionNums(k), sessionDir);
            if length(whichbetas)==1, fname = [fname,'_', condNames{whichbetas}]; end
            disp(fname);
            figurewrite(fname,[],[],'megfigs',1);
        elseif strcmp(plotType,'noise')
            figurewrite(sprintf('noisepool_%s',sessionDir),[],[],'megfigs',1);
        end
    end
end

%% Plot SNR before and after - all subjects 

axismin = 0; axismax = 20;
colors = jet(12);
allNoisepool = logical(allNoisepool);
hold on;
for nn = 1:length(sessionNums)
    fprintf(' session %d \n', sessionNums(nn));
    
    subplot(2,3,nn); cla;
    plot(allSNR1(nn,:),allSNR2(nn,:),'o','color',[1,1,1]*0.2); hold on;
    plot(allSNR1(nn,allNoisepool(nn,:)),allSNR2(nn,allNoisepool(nn,:)),'or');
    axismax = max([allSNR1(nn,:),allSNR2(nn,:)])*1.2;
    line([axismin,axismax],[axismin,axismax],'color','k');
    xlim([axismin,axismax]); ylim([axismin,axismax]); axis square;
    %xlabel('orig model SNR'); ylabel('final model SNR');
    title(sprintf('S%d', sessionNums(nn)))
    makeprettyaxes(gca,14,12);
end

if printFigsToFile
    figurewrite('SNRbeforeafter_allsubjs',[],[],'megfigs',1);
end

%% SNR before and after, mean across good PCs
% depends on having computed allPCchan
mean_before = zeros(1,length(sessionNums));
mean_after  = zeros(1,length(sessionNums));
for nn = 1:length(sessionNums)
    snr1 = allSNR1(nn,~isnan(allSNR1(nn,:)));
    snr2 = allSNR2(nn,~isnan(allSNR2(nn,:)));
    disp([length(snr1), length(snr2)]);
    mean_before(nn) = mean(snr1(logical(allPCchan{nn})),2);
    mean_after(nn)  = mean(snr2(logical(allPCchan{nn})),2);
end
colors = copper(5);
tmp = [mean_before; mean_after]';
hold on;
for nn = 1:length(sessionNums)
    plot(1:2,tmp(nn,:),'o-','color',colors(nn,:),'linewidth',2);
end

legend({'S1','S2','S3','S4','S5'});
xlim([0,3]);
set(gca,'xtick',1:2,'xticklabel',{'Before','After'}); 
ylabel('SNR');
makeprettyaxes(gca,14,14); axis square;

if printFigsToFile
    figurewrite('SNRbeforeafter_allsubjs2',[],[],'megfigs',1);
end