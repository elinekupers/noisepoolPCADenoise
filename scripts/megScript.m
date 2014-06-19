clear all;
inputDataDir = '/Volumes/HelenaBackup/denoisesuite/tmpmeg/';
outputFigDir = 'megfigs';
sessionNums  = 1:8;
sensorDataStr = 'b2';    % input data file string 
fitDataStr    = [sensorDataStr,'fr_fitfull75']; % fit data file string
whichfun     = 2;        % which fit (usually only 1)

%%
printFigsToFile = true;

% what to plot 
pp.plotPCselectByR2= true;  % R^2 as a function of number of PCs
    pp.PCSelMethod = 'snr'; 
pp.plotbbMap       = false;  % Broadband activity before and after denoising

pp.plotbbMap2      = false; % same as above but slightly different format (includes SL)
    pp.plotbbType  = 'SNR'; % specifies datatype for plotbbMap2 (options: 'S','N','SNR','R2')
    pp.plotbbConds = 1:3;   % conditions for plotbbMap2 (1=FULL,2=RIGHT,3=LEFT,1:3=All)
    
pp.plotNoisePool   = false; % location of noise pool 
pp.plotBeforeAfter = false;  % S, N, and SNR before and after denoising (all subjects togther)
    pp.doTop10     = true;  % specify format for plotBeforeAfter (top 10 or non-noise)
    
pp.plotSpectrum    = false; % spectrum of each channel, before and after denoising
pp.plotPCSpectrum  = false; % spectrum of PCs
pp.plotPCWeights   = false; % 

pp.condNames  = {'FULL','RIGHT','LEFT','OFF'};
pp.condColors = [0.1 0.1 0.9; 0.9 0.1 0.1; 0.1 0.9 0.1; .6 .6 .6];

%%
for k = 1:length(sessionNums)
    fprintf(' session %d \n', sessionNums(k));
    sessionDir = megGetDataPaths(sessionNums(k));
    % load fit file
    thisfile = fullfile(inputDataDir,sprintf('%s%s',sessionDir,fitDataStr));
    disp(thisfile); load(thisfile,'results');
    if pp.plotPCselectByR2, load(thisfile,'evalout'); end
    if pp.plotSpectrum, load(thisfile,'denoisedts'); end
    % load datafile
    datafile = fullfile(inputDataDir,sprintf('%s%s',sessionDir,sensorDataStr));
    disp(datafile); load(datafile,'badChannels');
    if pp.plotPCSpectrum || pp.plotSpectrum || pp.plotPCWeights, load(datafile,'sensorData','design'); end
    
    % some variables from the fit
    noisepool = results.noisepool;
    opt = results.opt;
    npcs2try = opt.npcs2try;
    %% ----------------------------------------------------
    %----------------------------------------------------
    % look at R2 as a function of PCs
    if pp.plotPCselectByR2
        
        %[chosen,pcchan,xvaltrend] = getpcchan(evalout(:,whichfun),noisepool,10,1.05,pp.PCSelMethod);
        switch pp.PCSelMethod
            case 'r2'
                metric = cat(1,evalout(:,whichfun).r2);
            case 'snr'
                metric = max(abs(cat(3,evalout(:,whichfun).beta_md)),[],1) ./ mean(cat(3,evalout(:,whichfun).beta_se),1);
                metric = squeeze(metric)';
        end
        chosen = results.pcnum(whichfun);
        xvaltrend = mean(metric(:,results.pcchan{whichfun}),2);
        
        % set up figure
        if k == 1, h1 = figure('position',[1,600,400,800]); end, figure(h1);
        % plot
        subplot(2,1,1);
        plot(0:npcs2try, metric);
        title(sprintf('%s for individual channels', upper(pp.PCSelMethod)))
        
        subplot(2,1,2);
        plot(0:npcs2try, xvaltrend(:,1), 'k','linewidth',2);
        title(sprintf('PC = %d', chosen(1)));
        
        for ii = 1:2
            subplot(2,1,ii);
            %xlabel('n pcs'); ylabel('R2');
            axis square;
            xlim([0,npcs2try]); %tmp = get(gca,'ylim'); ylim([-1,tmp(2)]);
            makeprettyaxes(gca,14);
        end
        vline(chosen,'k');
        hh = suptitle(sprintf('N%d : %s', sessionNums(k), sessionDir));
        set(hh,'interpreter','none');
        % write file
        if printFigsToFile
            figurewrite(sprintf('R2%02d_%s%s', sessionNums(k), sessionDir, fitDataStr),[],[], outputFigDir, 1);
        else
            pause;
        end
    end
    
    %% ----------------------------------------------------
    %----------------------------------------------------
    % look at broadband activity on the scalp
    if pp.plotbbMap
        % set up figure
        if k == 1, h2 = figure('position',[1,600,1200 400]); end, figure(h2);
        plotbbSNR(results, badChannels, pp.plotbbConds, whichfun, h2, pp.plotbbType);
        
        % write file
        if printFigsToFile
            figname = sprintf('%s%02d_%s%s', pp.plotbbType,sessionNums(k),sessionDir,fitDataStr);
            if length(pp.plotbbConds)==1, figname = [figname,'_', pp.condNames{pp.plotbbConds}]; end
            figurewrite(figname,[],[],outputFigDir,1);
        else
            pause;
        end
    end
    
    %% ----------------------------------------------------
    %----------------------------------------------------
    % look at broadband activity on the scalp, in comparison to stimulus
    % locked
    if pp.plotbbMap2
        % load stimulus locked
        slresults = load(fullfile(inputDataDir,sprintf('%s%sfSL_fitfull75',sessionDir,sensorDataStr)));
        slresults = slresults.results;
        % set up figure
        if k == 1, h3 = figure('position',[1,600,1400,400]); end, figure(h3);
        % plot
        switch pp.plotbbType
            case {'SNR','S','N'}
                sl_snr1 = getsignalnoise(slresults.origmodel(1),pp.plotbbConds, pp.plotbbType);
                ab_snr1 = getsignalnoise(results.origmodel(whichfun),  pp.plotbbConds, pp.plotbbType);
                ab_snr2 = getsignalnoise(results.finalmodel(whichfun), pp.plotbbConds, pp.plotbbType);
                clims_sl = [0, max(sl_snr1)];
                clims_ab = [0, max([ab_snr1, ab_snr2])];
            case 'R2'
                sl_snr1 = slresults.origmodel(1).r2;
                ab_snr1 = results.origmodel(whichfun).r2;
                ab_snr2 = results.finalmodel(whichfun).r2;
                clims_sl = [min(sl_snr1), max(sl_snr1)];
                clims_ab = [min([ab_snr1, ab_snr2]), max([ab_snr1, ab_snr2])];
        end
        sl_snr1a = to157chan(sl_snr1,~badChannels,'nans');
        ab_snr1a = to157chan(ab_snr1,~badChannels,'nans');
        ab_snr2a = to157chan(ab_snr2,~badChannels,'nans');
        
        subplot(1,3,1);
        megPlotMap(sl_snr1a,clims_sl,h3,'jet','Stimulus Locked Original');
        subplot(1,3,2);
        megPlotMap(ab_snr1a,clims_ab,h3,'jet','Broad Band Original');
        subplot(1,3,3);
        megPlotMap(ab_snr2a,clims_ab,h3,'jet',sprintf('Broad Band PC %d',results.pcnum(whichfun)));
        
        % write file
        if printFigsToFile
            figname = sprintf('%sMap_%02d_%s%s', pp.plotbbType,sessionNums(k),sessionDir,fitDataStr);
            if length(pp.plotbbConds)==1, figname = [figname,'_', pp.condNames{pp.plotbbConds}]; end
            figurewrite(figname,[],[],outputFigDir,1);
        else
            pause;
        end
    end
    
    %% ----------------------------------------------------
    %----------------------------------------------------
    % Look at the location of the noise pool
    if pp.plotNoisePool
        noisepool = results.noisepool;
        noise2   = to157chan(noisepool,~badChannels,0);
        % plot 
        figure(100); 
        megPlotMap(noise2,[0,1],[],'autumn',sprintf('Noise channels: N = %d',sum(noisepool)));
        colorbar off; 
        % write file 
        if printFigsToFile
            figurewrite(sprintf('noisepool_%s',sessionDir),[],[], outputFigDir, 1);
        else
            pause;
        end
    end
    
    %% ----------------------------------------------------
    %----------------------------------------------------
    % Look at SNR before and after 
    if pp.plotBeforeAfter
        if pp.doTop10
            pcchan = results.pcchan{whichfun}; st = 'Top 10';
        else
            pcchan = ~results.noisepool; st = 'Non Noise'; 
        end
        ab_signal1 = abs(results.origmodel(whichfun).beta_md(:,pcchan));
        ab_noise1  = results.origmodel(whichfun).beta_se(:,pcchan);
        ab_signal2 = abs(results.finalmodel(whichfun).beta_md(:,pcchan));
        ab_noise2  = results.finalmodel(whichfun).beta_se(:,pcchan);
        ab_snr1    = ab_signal1./ab_noise1;
        ab_snr2    = ab_signal2./ab_noise2;
         
        % set up figure
        c = ['b','r','g'];
        if k == 1, h4 = figure('position',[1,600,1200,600]); end, figure(h4);
        
        % plot 
        % signal
        subplot(3,length(sessionNums),k); cla; hold on; 
        for nn = 1:3
            plot(ab_signal1(nn,:),ab_signal2(nn,:),'o','color',pp.condColors(nn,:));
        end
        axis square;
        axismax = max([ab_signal1(:);ab_signal2(:)])*1.2;
        xlim([0,axismax]); ylim([0,axismax]); line([0,axismax],[0,axismax],'color','k');
        title(sprintf('S%d : signal', sessionNums(k)));
        % noise
        subplot(3,length(sessionNums),k+length(sessionNums)); cla; hold on;
        for nn = 1:3
            plot(ab_noise1(nn,:),ab_noise2(nn,:),'o','color',pp.condColors(nn,:));
        end
        axismax = max([ab_noise1(:); ab_noise2(:)])*1.2;
        xlim([0,axismax]); ylim([0,axismax]); line([0,axismax],[0,axismax],'color','k');
        axis square;
        title(sprintf('S%d : noise', sessionNums(k)));
        % snr
        subplot(3,length(sessionNums),k+2*length(sessionNums)); cla; hold on;
        for nn = 1:3
            plot(ab_snr1(nn,:),ab_snr2(nn,:),'o','color',pp.condColors(nn,:));
        end
        axismax = max([ab_snr1(:); ab_snr2(:)])*1.2;
        xlim([0,axismax]); ylim([0,axismax]); line([0,axismax],[0,axismax],'color','k');
        axis square;
        title(sprintf('S%d : SNR', sessionNums(k)));
        drawnow; sum(pcchan&(~results.noisepool))
        
        if k == length(sessionNums)
            for kk = 1:length(sessionNums)*3
                subplot(3,length(sessionNums),kk);
                %xlabel('orig model'); ylabel('final model');
                makeprettyaxes(gca,12);
            end
            suptitle(sprintf('Original versus Final Model : %s',st));
            if printFigsToFile
                figurewrite(sprintf('SignalNoise_allsubjs%s',fitDataStr),[],[], outputFigDir, 1);
            end
        end
    end
    
    
    %% ----------------------------------------------------
    %----------------------------------------------------
    % power spectrum 
    if pp.plotSpectrum
        
        if k == 1, h5 = figure('Position',[0,600,1200,500]); end, figure(h5);
        for chanNum = 1:157
            chanNum0 = megGetOrigChannel(chanNum,badChannels);
            if isnan(chanNum0), continue; end
            % results.origmodel.beta_md(:,chanNum0)
            % results.finalmodel.beta_md(:,chanNum0)
            for icond = 1:3
                epochConds = {design(:,icond)==1, all(design==0,2)};
                ax1 = subplot(1,2,1); cla;
                megPlotLogSpectra(sensorData,epochConds, badChannels, chanNum, ax1, pp.condNames([icond,4]));
                text(12,1e4,sprintf('beta=%0.2f',results.origmodel.beta_md(icond,chanNum0)),'fontsize',14,'color','r')
                
                ax2 = subplot(1,2,2); cla;
                megPlotLogSpectra(denoisedts{whichfun}, epochConds, badChannels, chanNum, ax2, pp.condNames([icond,4]));
                text(12,1e4,sprintf('beta=%0.2f',results.finalmodel.beta_md(icond,chanNum0)),'fontsize',14,'color','r')
                
                title(sprintf('PC = %d', results.pcnum(whichfun)));
                if printFigsToFile
                    figname = sprintf('spec%s_ch%03d_%s',sessionDir,chanNum,pp.condNames{icond});
                    figurewrite(figname,[],[],sprintf('%s/s%d',outputFigDir,sessionNums(k)),1);
                else
                    pause;
                end
            end
        end
    end
    
    
    %% ----------------------------------------------------
    %----------------------------------------------------
    % power spectrum of PCs
    if pp.plotPCSpectrum
        if ~isfield(results,'pcs')
            error('PCs not saved into results');
        end
        % get pcs into the right format 
        pcs = catcell(3,results.pcs); 
        pcs = permute(pcs,[2,1,3]);   % npcs x time x epochs
        epoch_idx = {design(:,1)==1, design(:,2)==1, design(:,3)==1, all(design==0,2)};
        
        % define figure properties
        if k == 1 
            f = (0:999);
            xl = [8 200];
            fok = f;
            fok(f<=xl(1) | f>=xl(2) ...
                | mod(f,60) < 1 | mod(f,60) > 59 ...
                ) = [];
            xt = [12:12:72, 96,144,192];
            h6 = figure('Position',[0,600,700,500]);
        end
        
        % plot 
        figure(h6);
        for p = 1:size(pcs,1)
            spec = abs(fft(squeeze(pcs(p,:,:))))/size(pcs,2)*2;
            clf; hold on;
            for ii = 1:length(epoch_idx)
                this_data = spec(:,epoch_idx{ii}).^2;
                plot(fok, nanmean(this_data(fok+1,:),2),  '-',  'Color', pp.condColors(ii,:), 'LineWidth', 2);
            end
            legend(pp.condNames,'location','best');
            set(gca, 'XLim', xl, 'XTick', xt, 'XScale', 'log', 'FontSize', 12);
            ss = 12; yl = get(gca, 'YLim'); for ii =ss:ss:180, plot([ii ii], yl, 'k--'); end
            title(sprintf('PC number %d', p));
            
            if printFigsToFile
                figurewrite(sprintf('%s_PC%02d',sessionDir,p),[],[],sprintf('%s/s%d',outputFigDir,sessionNums(k)),1);
            else
                pause;
            end
        end        
    end
    
    %% ----------------------------------------------------
    %----------------------------------------------------
    % Weights of PCs for regression 
    if pp.plotPCWeights
        
        if ~isfield(results,'pcs')
            error('PCs not saved into results');
        end
        
        % get pcs into the right format 
        pcs = catcell(3,results.pcs); 
        pcs = permute(pcs,[2,1,3]);   % npcs x time x epochs
        epoch_idx = {design(:,1)==1, design(:,2)==1, design(:,3)==1, all(design==0,2)};
        % get weights 
        fprintf('computing weights...\n');
        nepoch = size(sensorData,3);
        wts = [];
        for rp = 1:nepoch
            currsig = sensorData(:,:,rp)';
            wts = cat(3,wts,results.pcs{rp}(:,1:results.pcnum(whichfun))\currsig); % pcnum x channum
        end
        % average across epochs per condition
        fprintf('averaging weights...\n');
        wts_mean = cell(1,length(epoch_idx));
        for ii = 1:length(epoch_idx)
            wts_mean{ii} = mean(wts(:,:,epoch_idx{ii}),3);
        end
        
        % define figure properties
        if k == 1 
            h7 = figure('Position',[0,600,800,1000]);
            maxpcs = 20;
        end
        % plot weights for each condition 
        for ii = 1:4
            subplot(3,2,ii);
            if ii == 1
                imagesc(wts_mean{ii});
                clims = get(gca,'clim');
                clims = max(abs(clims))*[-1,1];
            end
            imagesc(wts_mean{ii}(1:maxpcs,:),clims); colorbar;
            title(pp.condNames{ii}); xlabel('Channel'); ylabel('PC number');
            makeprettyaxes(gca);
        end
        % average together as a function of PC number 
        subplot(3,2,[5,6]);
        hold on;
        for ii = 1:4
            plot(1:results.pcnum(whichfun), mean(wts_mean{ii}(:,~results.noisepool),2), ...
                'Color', pp.condColors(ii,:),'linewidth',2);
        end
        title('non noise channels'); xlabel('PC number'); ylabel('average weight');
        makeprettyaxes(gca);
        % print to file 
        if printFigsToFile
            figurewrite(sprintf('%02d_%s_PCWeight',k,sessionDir),[],[],sprintf('%s/s%d',outputFigDir,sessionNums(k)),1);
        else
            pause;
        end
    end 
end
