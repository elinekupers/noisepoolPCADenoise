function fH = plotSNRPrePostPanelABC(dataAll,exampleSessions,condColors,figureDir,saveFigures,figurenumber)

fH = figure; set(fH, 'Color', 'w');

for k = 1:numel(exampleSessions)

    dataAllind = ~cellfun(@isempty,dataAll);
    dataAll = dataAll(dataAllind);
    % top 10 channels
    pcchan = getTop10(dataAll{k}{1}.results);
    
    % signal and noise before denoising
    ab_signal1 = abs(dataAll{k}{1}.results.origmodel(1).beta_md(:,pcchan));
    ab_signal1 = dataAll{k}{1}.results.origmodel(1).beta_md(:,pcchan);

    ab_noise1  = dataAll{k}{1}.results.origmodel(1).beta_se(:,pcchan);
    % signal and noise after denoising 
    ab_signal2 = abs(dataAll{k}{1}.results.finalmodel(1).beta_md(:,pcchan));
    ab_signal2 = dataAll{k}{1}.results.finalmodel(1).beta_md(:,pcchan);

    ab_noise2  = dataAll{k}{1}.results.finalmodel(1).beta_se(:,pcchan);
    % snr = signal/denoise
    ab_snr1    = ab_signal1./ab_noise1;
    ab_snr2    = ab_signal2./ab_noise2;
    
    % plot each condition as a different color 
    % signal 
    subplot(3,length(exampleSessions),k); cla; hold on;
    for nn = 1:3
        plot(ab_signal1(nn,:),ab_signal2(nn,:),'o','color',condColors(nn,:));
    end
    axis square;
    axismax = max([ab_signal1(:);ab_signal2(:)])*1.2;
    xlim([0,axismax]); ylim([0,axismax]); line([0,axismax],[0,axismax],'color','k');
    title(sprintf('S%d : signal', exampleSessions(k)));
    makeprettyaxes(gca,9,9);
    
    % noise
    subplot(3,length(exampleSessions),k+length(exampleSessions)); cla; hold on;
    for nn = 1:3
        plot(ab_noise1(nn,:),ab_noise2(nn,:),'o','color',condColors(nn,:));
    end
    axismax = max([ab_noise1(:); ab_noise2(:)])*1.2;
    xlim([0,axismax]); ylim([0,axismax]); line([0,axismax],[0,axismax],'color','k');
    axis square;
    title(sprintf('S%d : noise', exampleSessions(k)));
    makeprettyaxes(gca,9,9);
    
    % snr
    subplot(3,length(exampleSessions),k+2*length(exampleSessions)); cla; hold on;
    for nn = 1:3
        plot(ab_snr1(nn,:),ab_snr2(nn,:),'o','color',condColors(nn,:));
    end
    axismax = max([ab_snr1(:); ab_snr2(:)])*1.2;
    xlim([0,axismax]); ylim([0,axismax]); line([0,axismax],[0,axismax],'color','k');
    axis square;
    title(sprintf('S%d : SNR', exampleSessions(k)));
    makeprettyaxes(gca,9,9);
    
    drawnow;
end

if saveFigures
    if isequal(figurenumber,7)
        figurewrite(fullfile(figureDir,'Figure7abc_snrexamples'),[],0,'.',1);
    elseif isequal(figurenumber,11)
        figurewrite(fullfile(figureDir,'Figure11abc_snrexamples'),[],0,'.',1);
    end
end