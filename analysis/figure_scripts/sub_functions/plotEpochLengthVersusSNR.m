function fH = plotEpochLengthVersusSNR(whichSubjects,dataAll,epochDurs,condColors,saveFigures,figureDir)

npsIdx = 2;
snr_diff = zeros(length(whichSubjects),length(epochDurs),3);
nepochs  = zeros(1,length(whichSubjects));

for whichSubject = whichSubjects

    
    % get all the results
    results_all = catcell(1,dataAll{whichSubject}{1}.allResults(:,npsIdx));
    
    % get top 10 channels and use the same 10 channels for all epoch
    % durations
    pcchan = getTop10(results_all(1));
    
    % total number of epochs for this session 
    nepochs(whichSubject) = length(results_all(1).opt.epochgroup);
    
    % compute the difference between pre and post
    for nn = 1:length(epochDurs)
        % get top 10 channels. note the choice of top 10 changes for each glm
        % result returned (for each epoch duration)
        %pcchan = getTop10(results_all(nn),whichfun);
        for icond = 1:3
            snr_pre  = getsignalnoise(results_all(nn).origmodel,icond);
            snr_post = getsignalnoise(results_all(nn).finalmodel,icond);
            snr_diff(whichSubject,nn,icond) = mean(snr_post(pcchan)-snr_pre(pcchan));
        end
    end
end

%%
% define colors - vary saturation across individual subjects
satValues = linspace(0.1,1,8);
colorRGB = varysat(condColors,satValues);

% set up figure and plot 
fH = figure('position',[0,300,250,600]);
for icond = 1:3 % for each condition 
    subplot(3,1,icond); hold on;
    % plot snr difference versus epoch duration 
    for nn = 1:length(whichSubjects)
        plot([epochDurs(1:end-1),nepochs(nn)], snr_diff(nn,:,icond),'o-','color',squeeze(colorRGB(icond,nn,:)),'linewidth',2);
    end
    % format axes and make figure pretty 
    set(gca,'xscale','log');
    xlim([0.5,1500]); set(gca,'xtick',epochDurs,'xscale','log');
    ylim([-3,7]); set(gca, 'YTick', [-2:2:7])
    ylabel('Difference in SNR (post-pre)');
    makeprettyaxes(gca,9,9);
end

if saveFigures
    figurewrite(fullfile(figureDir,'figure_epochdur'),[],0,'.',1);
end