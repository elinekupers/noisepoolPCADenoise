function fH = plotNPCvsNoisepoolSNR(whichSubjects,dataAll,npools,npcs,saveFigures,figureDir)

%% Plot difference in SNR (post-pre) as a function of number of channels in
%% noise pool and number of PCs removed

% nchan in noisepool x npcs removed x conds x nsessions
allSubjectsResults = [];

for whichSubject = whichSubjects
    fprintf(' compute SNR for subject %d \n', whichSubject);    
    
    % load results
    allvals = nan(length(npools),length(npcs),3);
    for np = 1:length(npools)
        for nc = 1:length(npcs)
            results = dataAll{whichSubject}{1}.allResults(np,nc);
            results = results{1};
            if isempty(results), continue; end
            
                % get top 10
                pcchan = getTop10(results);
          
                bb_signal_pre = results.origmodel(1).beta_md(:,pcchan);
                bb_noise_pre  = results.origmodel(1).beta_se(:,pcchan);
            
            	bb_signal_post = results.finalmodel(1).beta_md(:,pcchan);
                bb_noise_post  = results.finalmodel(1).beta_se(:,pcchan);
            
                bb_snr_pre     = bb_signal_pre./bb_noise_pre;
                bb_snr_post    = bb_signal_post./bb_noise_post;
            
                
                allvals(np,nc,:) = mean((bb_snr_post-bb_snr_pre),2);
                
        end
    end
    % concate across sessions
    allSubjectsResults = cat(4,allSubjectsResults,allvals);
end

%%
fH = figure('position',[0,300,450,900],'color','w');
clims = [[0,4];[0,3];[0,3]];
conditionNames = {'FULL','RIGHT','LEFT'};
for icond = 1:3
    subplot(3,1,icond);
    imagesc(1:length(npcs),1:length(npools),mean(allSubjectsResults(:,:,icond,:),4),clims(icond,:));
    
    set(gca,'ydir','normal');
    xlabel('Number of PCs removed');
    ylabel('Number of Channels in Noise pool');
    
    makeprettyaxes(gca,9,9);
    set(gca,'xtick',1:length(npcs),'ytick',1:length(npools),...
        'xticklabel',cellstr(num2str(npcs','%d')),'yticklabel',cellstr(num2str(npools','%d')));
    axis image; ch = colorbar; makeprettyaxes(ch,9,9);
    title(conditionNames{icond});
end

if saveFigures
    printnice(gcf, [1 300], fullfile(figureDir),'SF2GridSubjMean_NPCSvsNoisePool');
end