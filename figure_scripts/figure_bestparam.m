%% Define paths and data sets 
inputDataDir = '/Volumes/HelenaBackup/denoisesuite/tmpmeg/';
sessionNums  = [11:12,3:6,9:10];
fitDataStr    = 'b2fr_hpf2_fitall';
whichfun      = 1;
figuredir = 'manuscript_figs/figure_controls';

%% Plot difference in SNR (post-pre) as a function of number of channels in
%% noise pool and number of PCs removed 

% nchan in noisepool x npcs removed x conds x nsessions 
allvals2 = [];

for k = 1:length(sessionNums)
    fprintf(' session %d \n', sessionNums(k));
    sessionDir = megGetDataPaths(sessionNums(k));
    
    % load results 
    thisfile = fullfile(inputDataDir,sprintf('%s%s',sessionDir,fitDataStr));
    disp(thisfile); load(thisfile,'allresults','npools','npcs');
    
    allvals = nan(length(npools),length(npcs),3);
    for np = 1:length(npools)
        for nc = 1:length(npcs)
            results = allresults(np,nc);
            if isempty(results.finalmodel), continue; end
            
            % get top 10
            pcchan = getTop10(results,whichfun);
            
            % get snr before and after, nconds x nchannels
            ab_signal1 = abs(results.origmodel(whichfun).beta_md(:,pcchan));
            ab_noise1  = results.origmodel(whichfun).beta_se(:,pcchan);
            ab_signal2 = abs(results.finalmodel(whichfun).beta_md(:,pcchan));
            ab_noise2  = results.finalmodel(whichfun).beta_se(:,pcchan);
            ab_snr1    = ab_signal1./ab_noise1;
            ab_snr2    = ab_signal2./ab_noise2;
            % get snr difference for each condition (post-pre)
            for icond = 1:3
                allvals(np,nc,icond) = mean(ab_snr2(icond,:))-mean(ab_snr1(icond,:));
            end
        end
    end
    % concate across sessions 
    allvals2 = cat(4,allvals2,allvals);
end

%% 
fH = figure('position',[0,300,300,600]);
clims = [[0,4];[0,2];[0,2]];
for icond = 1:3
    subplot(3,1,icond);
    imagesc(1:length(npcs),1:length(npools),mean(allvals2(:,:,icond,:),4),clims(icond,:));
    
    set(gca,'ydir','normal');
    xlabel('Number of PCs removed'); 
    ylabel('Number of Channels in Noise pool');
    
    makeprettyaxes(gca,9,9); 
    set(gca,'xtick',1:length(npcs),'ytick',1:length(npools),...
    'xticklabel',cellstr(num2str(npcs','%d')),'yticklabel',cellstr(num2str(npools','%d')));
    axis image; ch = colorbar; makeprettyaxes(ch,9,9);
end
%figurewrite(fullfile(figuredir,'figure_grid_subjmean'),[],0,'.',1);
