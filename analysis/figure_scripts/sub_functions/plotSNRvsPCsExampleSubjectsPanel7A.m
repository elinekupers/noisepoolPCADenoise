function fH = plotSNRvsPCsExampleSubjectsPanel7A(dataAll,exampleSessions,condColors,axmax,figureDir,saveFigures, ylims, plotstr)

%% SNR increase for each MEG channel as a function of number of PCs removed for one example subject -
%% Fig. 7A

n = length(dataAll);

% Define colors
linecolors=copper(size(dataAll{n}.results.origmodel.beta,2));

fH = figure('position',[1,200,600,200]); set(fH, 'Color', 'w', 'Name', plotstr, 'NumberTitle', 'off')
for k = 1:length(exampleSessions)
    
    dataAllind = ~cellfun(@isempty,dataAll);
    dataAll = dataAll(dataAllind);
    
    % Define colors
    linecolors=copper(size(dataAll{k}.results.origmodel.beta,2));
    
    
    % get snr
    snr = (cat(3,dataAll{k}.evalout(:,1).beta_mn)) ./ cat(3,dataAll{k}.evalout(:,1).beta_sd);
    
    % plot for each condition
    for icond = 1:3 % Full, left, right
        subplot(length(exampleSessions),3,(k-1)*3+icond);
        hold on;
        this_snr = squeeze(snr(icond,:,:))';
        % plot each channel's snr as a function of number of pc's
        for ic = 1:size(this_snr,2)
            plot(0:axmax,this_snr(1:axmax+1,ic),'color',linecolors(ic,:));
        end
        % plot snr change for top10 channels
        xvaltrend = mean(this_snr(:,dataAll{k}.results.pcchan{1}),2);
        
        % plot snr change for noisepool
        noisepooltrend = squeeze(mean(snr(icond,dataAll{k}.results.noisepool,:),2));
        
        
        plot(0:axmax, xvaltrend(1:axmax+1,:), 'color', condColors(icond,:), 'linewidth',2);
        plot(0:axmax, noisepooltrend(1:axmax+1,:), 'color', condColors(icond,:), 'linewidth',5,'linestyle',':');
        
        % plot(axmax+1, xvaltrend(51,:), 'o', 'color', condColors(icond,:));
        axis square; xlim([0,axmax]);
        ylim([ylims(1),ylims(2)]); % if SL: ylim([0,50]) 
        set(gca,'ytick',[0:5:ylims(2)],'XTick',[0:5:axmax,axmax+2],'XTickLabel', [0:5:axmax  75]); 
        makeprettyaxes(gca,9,9);
    end
end
if saveFigures
    figurewrite(fullfile(figureDir,'FigureSNRvPCsData'),[],0,'.',1);
end