function fH = plotSNRvsPCsAllSubjectsPanel6B(dataAll,condColors,axmax,figureDir,saveFigures)

%% SNR increase as a function of number of PCs removed for all subjects -
%% Fig. 6B

% get the trend for the top 10 channels of all sessions
% files might take a while to load!
snr_top10 = [];
for k = 1:numel(dataAll)
    sessionDir = DFDgetdatapaths(sessionNums(k),conditionNumbers,inputDataDir);
    % load fit file
    thisfile = fullfile(inputDataDir,'savedProcData',sprintf('%s%s',sessionDir,fitDataStr));
    disp(thisfile); load(thisfile,'results','evalout');
    
    snr = abs(cat(3,evalout(:,whichFun).beta_md)) ./ cat(3,evalout(:,whichFun).beta_se);
    
    xvaltrend = [];
    for icond = 1:3
        this_snr = squeeze(snr(icond,:,1:11))';
        xvaltrend = cat(2, xvaltrend, mean(this_snr(:,results.pcchan{whichFun}),2));
    end
    snr_top10 = cat(3,snr_top10,xvaltrend);
end

%% Plot them together

% define colors - vary saturation for different subjects
satValues = linspace(0.1,1,8);
colorRGB = varysat(condColors,satValues);
ttls = {'FULL','RIGHT','LEFT'};

% plot for each condition
for icond = 1:3
    subplot(1,3,icond);hold on;
    for nn = 1:8 % for each subject
        plot(0:axmax, squeeze(snr_top10(:,icond,nn)), 'color', squeeze(colorRGB(icond,nn,:)));
    end
    axis square; xlim([0,axmax]);
    if plotBb
        ylim([0,12]); set(gca,'ytick',0:5:10);
    else
        ylim([0,40]); set(gca,'ytick',0:10:40);
    end
    title(ttls{icond});
    makeprettyaxes(gca,9,9);
end

if saveFigures
    if plotBb
        figurewrite(fullfile(figureDir,'Figure6BSNRvPCsAllSubjsBB'),[],0,'.',1);
    else
        figurewrite(fullfile(figureDir,'Figure11BSNRvPCsAllSubjsSL'),[],0,'.',1);
    end
end
