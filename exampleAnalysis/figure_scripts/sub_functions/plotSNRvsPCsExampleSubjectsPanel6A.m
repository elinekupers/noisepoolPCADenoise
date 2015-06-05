function fH = plotSNRvsPCsExampleSubjectsPanel6A(dataAll,exampleSessions,condColors,axmax,figureDir,saveFigures) %#ok<STOUT>

% Define colors
linecolors=copper(157);

fH = figure('position',[1,100,600,800]); set(fH, 'Color', 'w');
for k = 1:length(exampleSessions)
    % get snr
    snr = abs(cat(3,dataAll{k}{1}.evalout(:,1).beta_md)) ./ cat(3,dataAll{k}{1}.evalout(:,1).beta_se);
    
    % plot for each condition
    for icond = 1:3
        subplot(8,3,(k-1)*3+icond); hold on;
        this_snr = squeeze(snr(icond,:,:))';
        % plot each channel's snr as a function of number of pc's
        for ic = 1:size(this_snr,2)
            plot(0:axmax,this_snr(1:axmax+1,ic),'color',linecolors(ic,:));
        end
        % plot snr change for top10 channels
        xvaltrend = mean(this_snr(:,dataAll{k}{1}.results.pcchan{1}),2);
        plot(0:axmax, xvaltrend(1:axmax+1,:), 'color', condColors(icond,:), 'linewidth',2);
        %plot(axmax+1, xvaltrend(51,:), 'o', 'color', condColors(icond,:));
        axis square; xlim([0,axmax]);
        ylim([0,15]); % if SL: ylim([0,50])
        makeprettyaxes(gca,9,9);
    end
end
if saveFigures
   figurewrite(fullfile(figureDir,'Figure6SNRvPCsExampleSubjectsBB'),[],0,'.',1);
end