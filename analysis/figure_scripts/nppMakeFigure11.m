function nppMakeFigure11()
%% Function to reproduce Figure 11 with SNR difference before and after 
% denoising control analyses
%
% nppMakeFigure11()
% 
% Eline Kupers, Helena X. Wang, Kaoru Amano, Kendrick N. Kay, David J.
% Heeger, Jonathan Winawer. (2018) A non-invasive, quantitative study of
% broadband spectral responses in human visual cortex (PLOS ONE. VOLUME.
% ISSUE. DOI.)
%
% This figure will show a bargraph of broadband SNR difference before and
% after denoising for MEG Denoise and various control analyses.
%
% This function assumes that data is downloaded with the nppDownloaddata
% function. 



%% Choices to make:
whichSubjects        = 1:8;
dataDir              = fullfile(nppRootPath, 'analysis', 'data');   % Where to save data?
figureDir            = fullfile(nppRootPath, 'analysis', 'figures');% Where to save images?
saveFigures          = true;   % Save figures in the figure folder?
colors               = nppGetColors(3);
numOfControls        = 6;


%% Prepare data for figure
snr_diff = zeros(length(whichSubjects),numOfControls+1,3); % All controls, plus original result for all three conditions
for k = 1:length(whichSubjects)

    whichSubject = whichSubjects(k);
    fprintf(' Load subject %d \n', whichSubject);
    data = prepareData(dataDir,whichSubject,11);
    
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
fH = figure('position',[0,300,700,300],'Name', 'Figure 11', 'NumberTitle', 'off');
% define what the different conditions are 
types = {'noisepool-PCA',...                  1
    'Order shuffled', ...                   2
    'Random Amplitude', ...                 3
    'Phase-scrambled',...                   4
    'Replace PCs with random values', ...   5
    'All channels in noise pool', ...        6
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
    
    % BAR PLOT mean and sem across subjects 
%     mn  = mean(snr_diff2(:,:,icond));
%     sem = std(snr_diff2(:,:,icond))/sqrt(length(whichSubjects));
%     bar(1:nnull, mn,'EdgeColor','none','facecolor',colors(icond,:)); hold on
%     errorbar2(1:nnull,mn,sem,1,'-','color',colors(icond,:));
%     
% %     Plot individual subjects data points
%     for ii = 1:nnull
%         plot(ii*ones(size(snr_diff2(:,ii,icond))),snr_diff2(:,ii,icond),'o',...
%             'MarkerEdgeColor','w','MarkerFaceColor', 'k');
%     end
    
    % Or a BOXPLOT
     boxplot(snr_diff2(:,:,icond), 'Colors',colors(icond,:), 'BoxStyle','outline', 'Widths',0.2)
    
    % format figure and make things pretty 
    set(gca,'xlim',[0.2,nnull+0.8],'ylim',[-3,6]);
    makeprettyaxes(gca,9,9);
    set(gca,'XTickLabel',types(neworder));
    set(gca,'XTickLabelRotation',45);
%     set(get(gca,'XLabel'),'Rotation',45); 
    ylabel('Difference in SNR (post-pre)')
end

% Statistics
thisComparisonColumn = 1; % Which column do you want to compare the rest of your data with?

    for otherColumn = find([1:size(snr_diff2,2)]~=thisComparisonColumn)
        for icond = 1:3
            
            % Get p value by bootstrapping
            p = 2*(.5-abs(.5-mean(bootstrp(10000, @median, snr_diff2(:,thisComparisonColumn,icond) - snr_diff2(:,otherColumn,icond) )>0)));
            outBoot(otherColumn,icond) = p;
            
            % Traditional statistics
            [h,p] = ttest(snr_diff2(:,thisComparisonColumn,icond),snr_diff2(:,otherColumn,icond));            
            out(otherColumn,icond) = p;          
        end
    end

disp(out)
disp(outBoot)


if saveFigures
    figurewrite(fullfile(figureDir,'figure11_control'),[],0,'.',1);
end