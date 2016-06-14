function fH = plotSNRvsPCsAllSubjectsPanel7B(dataAll,condColors,axmax,figureDir,saveFigures)

%% SNR increase as a function of number of PCs removed for all subjects -
%% Fig. 7B

% get the trend for the top 10 channels of all sessions
snr_top10 = [];

dataAllind = ~cellfun(@isempty,dataAll);
dataAll = dataAll(dataAllind);

for k = 1:numel(dataAll)

    snr = (cat(3,dataAll{k}.evalout(:,1).beta_md)) ./ cat(3,dataAll{k}.evalout(:,1).beta_se);    
    xvaltrend = [];
    for icond = 1:3
        this_snr = squeeze(snr(icond,:,[1:11 76]))';
        xvaltrend = cat(2, xvaltrend, mean(this_snr(:,dataAll{k}.results.pcchan{1}),2));
    end
    snr_top10 = cat(3,snr_top10,xvaltrend);
end

%% Plot them together

% define colors - vary saturation for different subjects
satValues = 1-linspace(0.1,1,numel(dataAll));
colorRGB = varysat(condColors,satValues);
ttls = {'FULL','RIGHT','LEFT'};

fH = figure('position',[1,200,600,200]); set(fH, 'Color', 'w');
% plot for each condition
for icond = 1:3
    subplot(1,3,icond);hold on;
    for nn = 1:numel(dataAll) % for each subject
        plot(0:axmax, squeeze(snr_top10(1:11,icond,nn)), 'color', squeeze(colorRGB(icond,nn,:)));
        
        % Plot for 75 pcs regressed out
        plot(12, squeeze(snr_top10(12,icond,nn)), 'o', 'Color', squeeze(colorRGB(icond,nn,:)), 'LineWidth', 2, 'MarkerSize', 9);
    end
    axis square; xlim([0,axmax+3]);
    ylim([-1,12]); set(gca,'ytick',[0:5:10],'XTick',[0:5:10,12],'XTickLabel', [0:5:10 75]); % for SL: ylim([0,40]); set(gca,'ytick',0:10:40);
    title(ttls{icond});
    makeprettyaxes(gca,9,9);   
    xlabel('Number of PCs')
    ylabel('Median SNR (top 10 channels)')
end

if saveFigures
    figurewrite(fullfile(figureDir,'Figure7BSNRvPCsAllSubjsBB'),[],0,'.',1);
end



