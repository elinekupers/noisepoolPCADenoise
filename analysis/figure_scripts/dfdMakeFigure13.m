function dfdMakeFigure13()

%% Function to reproduce Figure 13 (Spatialmap) across all subjects of CiNet dataset
%
% dfdMakeFigure13()
%
% Eline Kupers, Helena X. Wang, Kaoru Amano, Kendrick N. Kay, David J.
% Heeger, Jonathan Winawer. (YEAR) Broadband spectral responses in visual
% cortex revealed by a new MEG denoising algorithm.
% (JOURNAL. VOLUME. ISSUE. DOI.)
%
% This figure will show an interpolated spatial map of the SNR values in
% each channel for the stimulus locked signal, broadband signals before
% using the denoising algorithm. The three separate conditions (Full,
% left, right hemifield stimulation are shown separately).
%
% This function assumes that data is downloaded with the DFDdownloaddata
% function.

%% Choices to make:
figureDir           = fullfile(dfdRootPath, 'analysis', 'figures'); % Where to save images?
dataDir             = fullfile(dfdRootPath, 'analysis', 'data');    % Where to save data?
saveFigures     	= true;     % Save figures in the figure folder?
threshold           = 0;
numCols             = 4;
meshData            = cell(1,numCols);

%% Get denoising results for raw data
whichSubjectsRaw    = 9:12;                % Raw
meshData([1 2 4])     = dfdMakeFigure13AcrossSubjects(whichSubjectsRaw,figureDir,dataDir,saveFigures,threshold);


%% Get denoising results for tsss data
whichSubjectsTSSS   = [14,16,18,20];       % TSSS
meshData([1 3 5])     = dfdMakeFigure13AcrossSubjects(whichSubjectsTSSS,figureDir,dataDir,saveFigures,threshold);

%% Make bar graph

% Get top ten across conditions for each subject
snr_mn = zeros(length(whichSubjectsRaw),size(meshData,2)-1,3); % All broadband columns
for whichSubject = 1:size(whichSubjectsRaw,2)
    
    
    for ii = [whichSubjectsRaw(whichSubject),whichSubjectsTSSS(whichSubject)]
        thisSubjectIdx = ii==[whichSubjectsRaw(whichSubject),whichSubjectsTSSS(whichSubject)];
        % Get noisepool info
        dd = load(sprintf(fullfile(dataDir, 's%02d_denoisedData_bb.mat'),ii));
        badChannels = dd.badChannels;
        noisePool(thisSubjectIdx,:) = to157chan(dd.results.noisepool,~badChannels,1); clear dd
    end
    noisePoolAcrossDenoisingConditions(whichSubject,:) = any(noisePool);
    
    % Extract data
    theseData = zeros(size(meshData,2)-1,3,size(noisePool,2));
    for thisColumn = 2:size(meshData,2)
        theseData(thisColumn-1,:,:) = meshData{thisColumn}(1:3,:,whichSubject);
    end
    
    % Get max snr across the three conditions
    thissnr = max(reshape(theseData,[],size(noisePool,2)));
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

%% Plot figure
fH = figure('position',[0,300,700,300]);
% define what the different conditions are
types = {'Raw','TSSS','MEG Denoise','MEG Denoise + TSSS'}; %
colors = dfdGetColors(3);

plotOrder = [1 3 2 4];
newOrder = types(plotOrder);

satValues = 1-linspace(0.1,1,4);
colorSaturated = varysat(colors,satValues);

nnull = length(types);
for icond = 1:3
    subplot(1,3,icond); hold on;
    % mean and sem across subjects
    %     mn  = mean(snr_mn(:,:,icond));
    %     sem = std(snr_mn(:,:,icond))/sqrt(4);
    %     bar(1:nnull, mn,'EdgeColor','none','facecolor',colors(icond,:)); hold on
    %     errorbar2(1:nnull,mn,sem,1,'-','color',colors(icond,:));
   for whichSubject = 1:4
        plot([1 2], snr_mn(whichSubject,[1 3],icond), 'o-', 'Color', squeeze(colorSaturated(icond,whichSubject,:)), ...
            'MarkerEdgeColor',squeeze(colorSaturated(icond,whichSubject,:)),'MarkerFaceColor', squeeze(colorSaturated(icond,whichSubject,:)), 'LineWidth', 2)
        
        plot([3 4], snr_mn(whichSubject,[2 4],icond), 'o-', 'Color', squeeze(colorSaturated(icond,whichSubject,:)), ...
            'MarkerEdgeColor',squeeze(colorSaturated(icond,whichSubject,:)),'MarkerFaceColor', squeeze(colorSaturated(icond,whichSubject,:)), 'LineWidth', 2)
   end
    % format figure and make things pretty
    set(gca,'xlim',[0.2,nnull+0.8],'ylim',[-3,12]);
    makeprettyaxes(gca,9,9);
    set(gca,'XTickLabel',newOrder, 'XTick', 1:4);
    set(gca,'XTickLabelRotation',45);
    %     set(get(gca,'XLabel'),'Rotation',45);
    ylabel('Broadband SNR')
end

if saveFigures
    figurewrite(fullfile(figureDir,'figure13b_linegraphTSSS'),[],0,'.',1);
end