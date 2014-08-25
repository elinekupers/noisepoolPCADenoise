%% Define paths and data sets 

inputDataDir = '/Volumes/HelenaBackup/denoisesuite/tmpmeg/';
fitDataStr   = 'b2fr_hpf2_fit10p1k_Nulls'; % fit data file string
ctrDataStr   = {'b2fr_hpf2_fit10Allnoise'};% additional controls to add 
sessionNums =  [11,12, 3:6, 9:10];%[1:6,9,10];
whichfun = 1;
condColors   = [63, 121, 204; 228, 65, 69; 116,183,74]/255;
figuredir = 'manuscript_figs/figure_controls';

%% Plot Different types of control analyses - Fig. 10 

snr_diff = zeros(length(sessionNums),5+length(ctrDataStr),3);
for k = 1:length(sessionNums)
    sessionDir = megGetDataPaths(sessionNums(k));
    % load fit file - may take a while 
    thisfile = fullfile(inputDataDir,sprintf('%s%s',sessionDir,fitDataStr));
    disp(thisfile); load(thisfile);
    
    % get all the results
    results_null = catcell(1,allResults);
    
    % load additional controls
    for jj = 1:length(ctrDataStr)
        thisfile = fullfile(inputDataDir,sprintf('%s%s',sessionDir,ctrDataStr{jj}));
        disp(thisfile); tmp = load(thisfile,'results');
        results_null = cat(1,results_null,tmp.results);
    end
    
    % get top 10 channels
    pcchan = getTop10(results_null(1),whichfun);
    
    % compute the difference between pre and post
    for icond = 1:3
        for nn = 1:length(results_null)
            snr_pre  = getsignalnoise(results_null(nn).origmodel(whichfun),icond);
            snr_post = getsignalnoise(results_null(nn).finalmodel(whichfun),icond);
            snr_diff(k,nn,icond) = mean(snr_post(pcchan)-snr_pre(pcchan));
        end
    end
end

%% 
fH = figure('position',[0,300,700,200]);
% define what the different conditions are 
types = {'Original','Phase-scrambled','Order shuffled', 'Random Amplitude', 'Random PCs','All Noise'};
% re-arrange the order of the bars 
neworder = [1,3,4,2,5,6];
newtypes = types(neworder);

snr_diff2 = snr_diff(:,neworder,:);
nnull = length(types);
% fH = figure('position',[0,300,200,200]);
for icond = 1:3
    subplot(1,3,icond);
    % mean and sem across subjects 
    mn  = mean(snr_diff2(:,:,icond));
    sem = std(snr_diff2(:,:,icond))/sqrt(8);
    bar(1:nnull, mn,'EdgeColor','none','facecolor',condColors(icond,:)); hold on
    errorbar2(1:nnull,mn,sem,1,'-','color',condColors(icond,:));
    % format figure and make things pretty 
    set(gca,'xlim',[0.2,nnull+0.8],'ylim',[-1,5]);
    makeprettyaxes(gca,9,9);
end

%figurewrite(fullfile(figuredir,'figure_allconds'),[],0,'.',1);