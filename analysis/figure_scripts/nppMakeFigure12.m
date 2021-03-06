function nppMakeFigure12
%% Function to reproduce Figure 12 (Spatialmap) across NYU dataset subjects having no, CALM or TSPCA preprocessing
%
% nppMakeFigure12()
%
% Eline Kupers, Helena X. Wang, Kaoru Amano, Kendrick N. Kay, David J.
% Heeger, Jonathan Winawer. (2018) A non-invasive, quantitative study of
% broadband spectral responses in human visual cortex (PLOS ONE. VOLUME.
% ISSUE. DOI.)
%
% This figure will show an interpolated spatial map of the SNR values in
% each channel for the stimulus locked signal, broadband signals before
% using the denoising algorithm. The three separate conditions (Full,
% left, right hemifield stimulation are shown separately).
%
% This function assumes that data is downloaded with the nppdownloaddata
% function.

%% Choices to make:
figureDir           = fullfile(nppRootPath, 'analysis', 'figures'); % Where to save images?
dataDir             = fullfile(nppRootPath, 'analysis', 'data');    % Where to save data?
saveFigures     	= true;     % Save figures in the figure folder?
threshold           = 0;
numCols             = 7;
meshData            = cell(1,numCols);
%% Get denoising results for raw (no preprocessing) data
% Makes columns 1,2,5
whichSubjectsRaw    = 1:8;                % Raw
meshData([1 2 5]) = nppMakeFigure12AcrossSubjects(whichSubjectsRaw,figureDir,dataDir,saveFigures,threshold);

%% Get denoising results for CALM preprocessed data
% Makes columns 3,6
whichSubjectsCALM   = 21:28;              % CALM
meshData([1 3 6]) = nppMakeFigure12AcrossSubjects(whichSubjectsCALM,figureDir,dataDir,saveFigures,threshold);

%% Get denoising results for TSPCA preprocessed data
% Makes columns 4,7
whichSubjectsTSPCA  = 29:36;              % TSPCA
meshData([1 4 7]) =  nppMakeFigure12AcrossSubjects(whichSubjectsTSPCA,figureDir,dataDir,saveFigures,threshold);

%% Make bar graph

% Get top ten across conditions for each subject
snr_mn = zeros(length(whichSubjectsRaw),size(meshData,2)-1,3); % All broadband columns
for whichSubject = 1:size(whichSubjectsRaw,2)
    
    
   for ii = [whichSubjectsRaw(whichSubject),whichSubjectsCALM(whichSubject),whichSubjectsTSPCA(whichSubject)]
       thisSubjectIdx = find(ii==[whichSubjectsRaw(whichSubject),whichSubjectsCALM(whichSubject),whichSubjectsTSPCA(whichSubject)]);
        
       % Get noisepool info
        dd = load(sprintf(fullfile(dataDir, 's%02d_denoisedData_bb.mat'),ii));
        badChannels = dd.badChannels; 
        noisePool(thisSubjectIdx,:) = to157chan(dd.results.noisepool,~badChannels,1); clear dd   
   end
   noisePoolAcrossDenoisingConditions(whichSubject,:) = any(noisePool);
    
   
   % Extract data
   theseData = zeros(size(meshData,2)-1,3,157);
   for thisColumn = 2:size(meshData,2)
       theseData(thisColumn-1,:,:) = meshData{thisColumn}(1:3,:,whichSubject);
   end

   % Get max snr across the three conditions
   thissnr = max(reshape(theseData,18,157));
   thissnr(noisePoolAcrossDenoisingConditions(whichSubject,:)) = -inf;
   [~,idx] = sort(thissnr,'descend');
   thispcchan = false(size(noisePoolAcrossDenoisingConditions(whichSubject,:)));
   thispcchan(idx(1:10))= 1;
    
   finalpcchan(whichSubject,:) = thispcchan;
    
    
    for thisColumn = 2:size(meshData,2)
        % compute the difference between pre and post
        for icond = 1:3
            snr_post  = meshData{thisColumn}(icond,finalpcchan(whichSubject,:),whichSubject);            
            snr_mn(whichSubject,thisColumn-1,icond) = mean(snr_post);
            
        end
    end
end

%% Statistics

thisComparisonColumn = 4;
for otherColumn = find([1:size(meshData,2)-1]~=thisComparisonColumn)
    for icond = 1:3
        % Get p value by bootstrapping
        p = 2*(.5-abs(.5-mean(bootstrp(10000, @median, snr_mn(:,thisComparisonColumn,icond) - snr_mn(:,otherColumn,icond) )>0)));
        outBoot(otherColumn,icond) = p;
        
%         Traditional statistics
        [h,p] = ttest(snr_mn(:,thisComparisonColumn,icond),snr_mn(:,otherColumn,icond));
        out(otherColumn,icond) = p;
    end
end
disp(outBoot)
disp(out)

%% Plot figure
fH = figure('position',[0,300,700,300],'Name', 'Figure 12B', 'NumberTitle', 'off');
% define what the different conditions are
types = {'Raw','CALM','TSPCA','noisepool-PCA','noisepool-PCA + CALM','noisepool-PCA + TSPCA'}; %
% re-arrange the order of the bars
colors = nppGetColors(3);

nnull = length(types);
for icond = 1:3
    subplot(1,3,icond);
    
    % BAR PLOT
    % mean and sem across subjects
%     mn  = mean(snr_mn(:,:,icond));
%     sem = std(snr_mn(:,:,icond))/sqrt(8);
%     bar(1:nnull, mn,'EdgeColor','none','facecolor',colors(icond,:)); hold on
%     errorbar2(1:nnull,mn,sem,1,'-','color',colors(icond,:));
    
    % Plot individual data points
%     for ii = 1:nnull
%         plot(ii*ones(size(snr_mn(:,ii,icond))),snr_mn(:,ii,icond),'o',...
%             'MarkerEdgeColor','w','MarkerFaceColor', 'k');
%     end
    
    % BOXPLOT
    boxplot(snr_mn(:,:,icond), 'Colors',colors(icond,:), 'BoxStyle','outline', 'Widths',0.2,'MedianStyle','target')
    
    
    % format figure and make things pretty
    set(gca,'xlim',[0.2,nnull+0.8],'ylim',[-1,10]);
    makeprettyaxes(gca,9,9);
    set(gca,'XTickLabel',types);
    set(gca,'XTickLabelRotation',45);
    %     set(get(gca,'XLabel'),'Rotation',45);
    ylabel('Broadband SNR')
end

if saveFigures
    % Only use figure write when producing high quality manuscript figure
    % (since it is very slow)
    hgexport(fH, fullfile(figureDir,'figure12b_bargraphCALMTSPCA'));
%     figurewrite(fullfile(figureDir,'figure12b_bargraphCALMTSPCA'),[],0,'.',1);
end