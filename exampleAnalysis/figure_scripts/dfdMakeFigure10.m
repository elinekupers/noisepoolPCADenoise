function dfdMakeFigure10()
%% Function to reproduce Figure 10 with SNR difference before and after 
% denoising control analyses
%
% dfdMakeFigure10()
% 
% AUTHORS. TITLE. JOURNAL. YEAR.
%
% This figure will show ...
%
% This function assumes that data is downloaded with the DFDdownloaddata
% function. 



%% Choices to make:
whichSubjects        = 1:8;
dataDir              = fullfile(dfdRootPath, 'exampleAnalysis', 'data');   % Where to save data?
figureDir            = fullfile(dfdRootPath, 'exampleAnalysis', 'figures');% Where to save images?
saveFigures          = true;   % Save figures in the figure folder?
condColors           = [63, 121, 204; 228, 65, 69; 116,183,74]/255;
num_of_controls      = 5;
dataAll              = [];

%% Load data for all subjects
for whichSubject = whichSubjects
    fprintf(' Load subject %d \n', whichSubject);
    [data,design,exampleIndex] = prepareData(dataDir,whichSubject,10);
    dataAll{whichSubject} = {data,design,exampleIndex}; %#ok<AGROW>
end 


%% Prepare data for figure
snr_diff = zeros(length(whichSubjects),num_of_controls+1,3); % All controls, plus original result for all three conditions
for k = 1:length(whichSubjects)
    
    results_null = [dataAll{k}{1}(1),dataAll{k}{1}{2}];
    
    % get top 10 channels 
    pcchan = getTop10(results_null{1,1}.results);
    
    % compute the difference between pre and post
    for icond = 1:3
        for nn = 1:length(results_null)
            snr_pre  = getsignalnoise(results_null{nn}.results.origmodel,icond);
            snr_post = getsignalnoise(results_null{nn}.results.finalmodel,icond);
            snr_diff(k,nn,icond) = mean(snr_post(pcchan)-snr_pre(pcchan));
        end
    end
end

%% Plot figure
fH = figure('position',[0,300,700,200]);
% define what the different conditions are 
types = {'Original','Phase-scrambled','Order shuffled', 'Random Amplitude', 'Random PCs','All Noise'};
% re-arrange the order of the bars 
neworder = [1,3,4,2,5,6];
newtypes = types(neworder);

snr_diff2 = snr_diff(:,neworder,:);
nnull = length(types);
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

if saveFigures
    figurewrite(fullfile(figureDir,'figure10_allconds'),[],0,'.',1);
end