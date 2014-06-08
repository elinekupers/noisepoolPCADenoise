clear all;
sessionNums = 3;
tmpmegdir = '/Volumes/HelenaBackup/denoisesuite/tmpmeg/';

%% plot R^2 as a function of the number of PCs
% load saved data 
% beware of long loading time

printFigsToFile = false;
optpcs = zeros(1,length(sessionNums));
plotType = 3;
loadNull = false;

if plotType == 1 % deprecated - compares mean of top10 to nonnoise/all
    figure('position',[1,600,800,300]);
elseif plotType == 2 % look at top 10 
    figure('position',[1,600,400,800]);
elseif plotType == 3 % look at spatial map of SNR 
    figure('position',[1,600,1200 400]);
end

ppfit = 'f_hpf2_fitfull75';
pp    = 'b2';
for k = 1:length(sessionNums)
    fprintf(' session %d \n', sessionNums(k));
    sessionDir = megGetDataPaths(sessionNums(k));
    % load fit file 
    thisfile = fullfile(tmpmegdir,sprintf('%s%s%s',sessionDir,pp,ppfit));
    %thisfile = fullfile(tmpmegdir,sprintf('%s_fitperm',sessionDir));
    disp(thisfile); load(thisfile,'results','evalout');
    % load data file 
    if plotType == 3
        load(fullfile(tmpmegdir,sprintf('%s%s',sessionDir,pp)),'badChannels'); 
    end

    fprintf(' done loading\n');
    if exist('results','var')
        noisepool = results.noisepool;
        opt = results.opt;
    end
    
    if plotType ~= 3
        % look at r2 as a function of pcs
        %-------------------------------------------
        [chosen,pcchan,xvaltrend] = getpcchan(evalout(:,1),noisepool,10,1.05);
        r2 = cat(1,evalout(:,1).r2);
        
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
        disp(opt.npcs); opt.npcs = size(xvaltrend,1)-1;
    end
    
    % plot
    if plotType == 1 % deprecated 
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
            xlim([0,opt.npcs]);
            vline(results.pcnum(1),'r'); hold on;
            %vline(finalmodel(1).pcnum,'r');
            vline(chosen,'k');
            makeprettyaxes(gca,12);
        end
        suptitle(sprintf('MEG Session %d : %s', sessionNums(k), sessionDir));
        
    elseif plotType == 2 % all channels versus top 10 

        subplot(2,1,1);
        plot(0:opt.npcs, r2(:,:,1));
        title('R^2 for individual channels')
        
        subplot(2,1,2);
        plot(0:opt.npcs, xvaltrend(:,1), 'k','linewidth',2);
        title(sprintf('PC = %d', chosen(1)));
        
        for ii = 1:2
            subplot(2,1,ii);
            %xlabel('n pcs'); ylabel('R2'); 
            axis square; 
            xlim([0,opt.npcs]); %tmp = get(gca,'ylim'); ylim([-1,tmp(2)]);
            makeprettyaxes(gca,14);
        end
        vline(chosen,'k');
        h2 = suptitle(sprintf('N%d : %s', sessionNums(k), sessionDir));
        set(h2,'interpreter','none');
    
    elseif plotType == 3 % SNR map 
        plotbbSNR(results,badChannels,1:3,1,gcf,'SNR');
    end
    
    %pause;
    %clear allEval results noisepool opt evalout
    fprintf('====================\n\n');
    
    if printFigsToFile
        figurewrite(sprintf('%02d_%s%s%s', sessionNums(k), sessionDir, pp, ppfit),[],[],'megfigs',1);
    end
end

%% Plot Comparisons with NULL
% plot R^2 as a function of PCs, with all subjects together
% also plot the different types of nulls 
figure('position',[1,600,1000,400]); 
fudge = [0,3,3,3,3];
ttls = {'Original','Phase scrambled','Order shuffled','Amplitude scrambled','Random pcs'};
plotType = 3;

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
        
    case 3
        ppfit = 'b2_epochGroup6s_fitfull75';
        for nn = 1:length(sessionNums)
            [sessionDir,megDataDir] = megGetDataPaths(sessionNums(nn));
            thisfile = fullfile(tmpmegdir,sprintf('%s%s_allnulls',sessionDir,ppfit));
            disp(thisfile); loaded = load(thisfile);
            
            r2 = loaded.allR2;
            thesechan = ~loaded.allResults{1}.noisepool;
            subplot(2,3,nn); cla; hold on;
            for k = 1:4
                plot(0:size(r2,1)-1, mean(r2(:,thesechan,k),2), colors(k));
            end
            
            axis square; xlim([0,75]); tmp = get(gca,'ylim'); ylim([-1,tmp(2)]);
            title(sprintf('S%d', sessionNums(nn)));
            makeprettyaxes(gca,14,12);
            if nn == length(sessionNums), legend(ttls,'location','bestoutside'); end 
        end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot the spatial maps as a function of denoising
fs = 16;
printFigsToFile = false;
optpcs = zeros(1,length(sessionNums));
plotType = 'noise';
plotSNR  = true;

if strcmp(plotType, 'SNR'), figure('position',[1,600,1600 500]);  end
allSNR1 = []; allSNR2 = []; allNoisepool = [];
whichbetas = 1:3; %<--- toggle here 

condNames = {'FULL','LEFT','RIGHT'};
pp = 'b2'; ppfit = 'f_hpf2_fitfull0';
for k = 1:length(sessionNums)
    fprintf(' session %d \n', sessionNums(k));
    [sessionDir,megDataDir,conditionNames] = megGetDataPaths(sessionNums(k), 1:6);
    %thisfile = fullfile('megfigs/matfiles',sprintf('%02d_%s',sessionNums(k),sessionDir));
    thisfile = fullfile(tmpmegdir,sprintf('%s%s%s',sessionDir,pp,ppfit));
    disp(thisfile); load(thisfile,'results');
    % load stimulus locked
    slresults = load(fullfile(tmpmegdir,sprintf('%s%s_stimlocked',sessionDir,pp)));
    slresults = slresults.results;
    % load bad channels
    load(fullfile(tmpmegdir,sprintf('%s%s',sessionDir,pp)),'badChannels'); 
    
    opt = results.opt; noisepool = results.noisepool;
    disp(opt.npoolmethod);
    
    if plotSNR
%         sl_snr1 = max(abs(results.origmodel(2).beta_md(whichbetas,:)),[],1)...
%             ./mean(results.origmodel(2).beta_se(whichbetas,:),1);
        sl_snr1 = max(abs(slresults.origmodel(1).beta_md(whichbetas,:)),[],1)...
            ./mean(slresults.origmodel(1).beta_se(whichbetas,:),1);
        ab_snr1 = max(abs(results.origmodel(1).beta_md(whichbetas,:)),[],1)./...
            mean(results.origmodel(1).beta_se(whichbetas,:),1);
        ab_snr2 = max(abs(results.finalmodel(1).beta_md(whichbetas,:)),[],1)./...
            mean(results.finalmodel(1).beta_se(whichbetas,:),1);
        clims_sl = [0, max(sl_snr1)];
        clims_ab = [0, max([ab_snr1, ab_snr2])];
        
    else
        sl_snr1 = slresults.origmodel(1).r2;
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot changes in Signal and Noise separately 

ppfit = 'b2f_hpf2_fitfull0';
condNames = {'FULL','LEFT','RIGHT'};
c = ['b','r','g']; 
figure('position',[1,600,1200,500]);
doTop10 = true;
for k = 1:length(sessionNums)
    fprintf(' session %d \n', sessionNums(k));
    [sessionDir] = megGetDataPaths(sessionNums(k));
    %thisfile = fullfile('megfigs/matfiles',sprintf('%02d_%s',sessionNums(k),sessionDir));
    thisfile = fullfile(tmpmegdir,sprintf('%s%s',sessionDir,ppfit));
    disp(thisfile); load(thisfile,'results');

    %pcchan = allPCchan{sessionNums(k)};
    if doTop10
        pcchan = results.opt.pcchan{1};
        st = 'Top 10';
    else
        pcchan = ~results.noisepool;
        st = 'Non Noise';
    end
    
    ab_signal1 = abs(results.origmodel(1).beta_md(:,pcchan));
    ab_noise1  = results.origmodel(1).beta_se(:,pcchan);
    ab_signal2 = abs(results.finalmodel(1).beta_md(:,pcchan));
    ab_noise2  = results.finalmodel(1).beta_se(:,pcchan);
    
    subplot(2,length(sessionNums),k); cla; hold on;
    for nn = 1:3
        plot(ab_signal1(nn,:),ab_signal2(nn,:),['o',c(nn)])
    end
    axis square;
    axismax = max([ab_signal1(:);ab_signal2(:)])*1.2;
    xlim([0,axismax]); ylim([0,axismax]); line([0,axismax],[0,axismax],'color','k');
    title(sprintf('S%d : signal', sessionNums(k)));
    xlabel('orig model'); ylabel('final model');
    makeprettyaxes(gca);
    
    subplot(2,length(sessionNums),k+length(sessionNums)); cla; hold on;
    for nn = 1:3
        plot(ab_noise1(nn,:),ab_noise2(nn,:),['o',c(nn)]);
    end
    axismax = max([ab_noise1(:); ab_noise2(:)])*1.2;
    xlim([0,axismax]); ylim([0,axismax]); line([0,axismax],[0,axismax],'color','k');
    axis square;
    title(sprintf('S%d : noise', sessionNums(k)));
    xlabel('orig model'); ylabel('final model');
    makeprettyaxes(gca);
    drawnow;
    %pause;
    sum(pcchan&(~results.noisepool))
end

h2=suptitle(st);
set(h2,'interpreter','none');

if printFigsToFile
    figurewrite(sprintf('SignalNoise_allsubjs%s',ppfit),[],[],'megfigs',1);
end

%% Look at the power spectrum of the PCs

k = 5; % session number 
% load data with pcs saved and corresponding design matrix 
[sessionDir,megDataDir] = megGetDataPaths(sessionNums(k));
thisfile = fullfile(tmpmegdir,sprintf('%sb2',sessionDir));
disp(thisfile); load(thisfile,'design','sensorData');
fitfile = fullfile(tmpmegdir,sprintf('%sb2_fitfull0',sessionDir));
disp(fitfile); load(fitfile,'results');
% get pcs into the right format 
pcs = catcell(3,results.pcs); 
pcs = permute(pcs,[2,1,3]);   % npcs x time x epochs
epoch_idx = {design(:,1)==1, design(:,2)==1, design(:,3)==1, all(design==0,2)};

%%
% define figure properties 
f = (0:999);
xl = [8 200];
fok = f; 
fok(f<=xl(1) | f>=xl(2) ...
    | mod(f,60) < 1 | mod(f,60) > 59 ...
    ) = [];
colors = [0.1 0.1 0.9; 0.9 0.1 0.1; 0.1 0.9 0.1; .6 .6 .6];
xt = [12:12:72, 96,144,192];
fH = figure('Position',[0,600,700,500]); 
for p = 1:size(pcs,1)
    spec = abs(fft(squeeze(pcs(p,:,:))))/size(pcs,2)*2;
    clf; hold on;
    for ii = 1:length(epoch_idx)
        this_data = spec(:,epoch_idx{ii}).^2;
        plot(fok, nanmean(this_data(fok+1,:),2),  '-',  'Color', colors(ii,:), 'LineWidth', 2);
    end
    legend('FULL','RIGHT','LEFT','OFF','location','best');
    set(gca, 'XLim', xl, 'XTick', xt, 'XScale', 'log', 'FontSize', 12);
    ss = 12; yl = get(gca, 'YLim'); for ii =ss:ss:180, plot([ii ii], yl, 'k--'); end
    title(sprintf('PC number %d', p));
    if printFigsToFile
        figurewrite(sprintf('%s_PC%02d',sessionDir,p),[],[],sprintf('megfigs/s%d',k),1);
    else
        pause;
    end
end

%% print PC weight separately per condition 
nepoch = size(sensorData,3);
wts = [];
for rp = 1:nepoch
    currsig = sensorData(:,:,rp)';
    wts = cat(3,wts,results.pcs{rp}(:,1:results.pcnum(1))\currsig); % pcnum x channum
end

figure('Position',[0,600,800,1000]);
ts = {'FULL','RIGHT','LEFT','OFF'};
for ii = 1:length(epoch_idx)
    wts_mean{ii} = mean(wts(:,:,epoch_idx{ii}),3);
end
maxpcs = 20;
for ii = 1:4
    subplot(3,2,ii); 
    if ii == 1
        imagesc(wts_mean{ii});
        clims = get(gca,'clim');
        clims = max(abs(clims))*[-1,1];
    end
    imagesc(wts_mean{ii}(1:maxpcs,:),clims);
    colorbar;
    title(ts{ii}); xlabel('Channel'); ylabel('PC number');
    makeprettyaxes(gca);
end
subplot(3,2,[5,6]);
hold on;
for ii = 1:4
    plot(1:results.pcnum(1), mean(wts_mean{ii}(:,~results.noisepool),2), ...
        'Color', colors(ii,:),'linewidth',2);
end
title('non noise channels'); xlabel('PC number'); ylabel('average weight');
makeprettyaxes(gca);
if printFigsToFile
    figurewrite(sprintf('%02d_%s_PCWeight',k,sessionDir),[],[],sprintf('megfigs/s%d',k),1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Look at power spectrum of a particular channel
for k = 1:6 % session number 
[sessionDir,megDataDir] = megGetDataPaths(sessionNums(k));
thisfile = fullfile(tmpmegdir,sprintf('%sb2',sessionDir));
disp(thisfile); load(thisfile,'sensorData','badChannels','design');
fitfile = fullfile(tmpmegdir,sprintf('%sb2f_epochGroup6so_fitfull75',sessionDir));
disp(fitfile); load(fitfile,'results','denoisedts');

load(fullfile(tmpmegdir,'epochGroups',[sessionDir, 'b2_epochGroup6so']));
discardepochs = epochGroup == 0;
sensorData = sensorData(:,:,~discardepochs);
design     = design(~discardepochs,:);

%%
fH = figure('Position',[0,600,1200,500]);
condNames = {'FULL','RIGHT','LEFT','OFF'};
for chanNum = 1:157
    chanNum0 = megGetOrigChannel(chanNum,badChannels);
    if isempty(chanNum0), continue; end
    % results.origmodel.beta_md(:,chanNum0)
    % results.finalmodel.beta_md(:,chanNum0)
    for icond = 1:3
        epochConds = {design(:,icond)==1, all(design==0,2)};
        
        ax1 = subplot(1,2,1); cla;
        megPlotLogSpectra(sensorData,epochConds, badChannels, chanNum, ax1, condNames([icond,4]));
        ax2 = subplot(1,2,2); cla;
        megPlotLogSpectra(denoisedts{1},epochConds, badChannels, chanNum, ax2,condNames([icond,4]));
        title(sprintf('PC = %d', results.pcnum(1)));
        if printFigsToFile
            figurewrite(sprintf('spec%s_ch%03d_%s',sessionDir,chanNum,condNames{icond}),[],[],sprintf('megfigs/s%d',k),1);
        else
            pause;
        end
    end
end

end