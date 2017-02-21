function dfdMakeFigure10()
%% Function to reproduce Figure 11 with SNR difference before and after 
% denoising control analyses
%
% dfdMakeFigure10()
% 
% Eline Kupers, Helena X. Wang, Kaoru Amano, Kendrick N. Kay, David J.
% Heeger, Jonathan Winawer. (YEAR) Broadband spectral responses in visual
% cortex revealed by a new MEG denoising algorithm.
% (JOURNAL. VOLUME. ISSUE. DOI.)
%
% This figure will show a bargraph of broadband SNR difference before and
% after denoising for MEG Denoise and various control analyses.
%
% This function assumes that data is downloaded with the DFDdownloaddata
% function. 



%% Choices to make:
whichSubjects        = 1:8;
dataDir              = fullfile(dfdRootPath, 'analysis', 'data');   % Where to save data?
figureDir            = fullfile(dfdRootPath, 'analysis', 'figures');% Where to save images?
saveFigures          = true;   % Save figures in the figure folder?
colors               = dfdGetColors(3);
numOfControls        = 6;


%% Prepare data for figure
snr_diff = zeros(length(whichSubjects),numOfControls+1,3); % All controls, plus original result for all three conditions
for k = 1:length(whichSubjects)

    whichSubject = whichSubjects(k);
    fprintf(' Load subject %d \n', whichSubject);
    data = prepareData(dataDir,whichSubject,10);
    
    results_null = [data(1),data{2}];
    
    % get top 10 channels 
    pcchan = getTop10(results_null{1,1}.results);
    
    % compute the difference between pre and post
    for icond = 1:3
        for nn = 1:length(results_null)
            the_contrast = [0 0 0];
            the_contrast(icond) = 1;
            snr_pre  = getsignalnoise(results_null{nn}.results.origmodel,the_contrast);
            snr_post = getsignalnoise(results_null{nn}.results.finalmodel,the_contrast);
            snr_diff(k,nn,icond) = mean(snr_post(pcchan)-snr_pre(pcchan));
        end
    end
end

%% Plot figure
fH = figure('position',[0,300,700,300]);
% define what the different conditions are 
types = {'MEG Denoise',...                  1
    'Order shuffled', ...                   2
    'Random Amplitude', ...                 3
    'Phase-scrambled',...                   4
    'Replace PCs with random values', ...   5
    'All channels in noisepool', ...        6
    'Concatenate epochs for denoising' ...  7
     }; % 
% re-arrange the order of the bars 
%neworder = [1,5,6];
neworder = [1, 7, 6, 4];
newtypes = types(neworder);

snr_diff2 = snr_diff(:,neworder,:);
nnull = length(neworder);
for icond = 1:3
    subplot(1,3,icond);
    % mean and sem across subjects 
    mn  = mean(snr_diff2(:,:,icond));
    sem = std(snr_diff2(:,:,icond))/sqrt(length(whichSubjects));
    bar(1:nnull, mn,'EdgeColor','none','facecolor',colors(icond,:)); hold on
    errorbar2(1:nnull,mn,sem,1,'-','color',colors(icond,:));
    % format figure and make things pretty 
    set(gca,'xlim',[0.2,nnull+0.8],'ylim',[-1,5]);
    makeprettyaxes(gca,9,9);
    set(gca,'XTickLabel',types(neworder));
    set(gca,'XTickLabelRotation',45);
%     set(get(gca,'XLabel'),'Rotation',45); 
    ylabel('Difference in SNR (post-pre)')
end

% Statistics
for thisColumn = 1:size(snr_diff2,2)
    for otherColumn = find([1:size(snr_diff2,2)]~=thisColumn)
        for icond = 1:3
            
            [h,p] = ttest(snr_diff2(:,thisColumn,icond),snr_diff2(:,otherColumn,icond));
            
            out(thisColumn,otherColumn,icond) = p;
            
        end
    end
end
disp(out)

if saveFigures
    figurewrite(fullfile(figureDir,'figure10_control'),[],0,'.',1);
end