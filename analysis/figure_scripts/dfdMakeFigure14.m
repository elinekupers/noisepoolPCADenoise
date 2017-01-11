function dfdMakeFigure14()

%% Function to reproduce Figure 14 (Spatialmap) across all subjects of CiNet dataset
%
% dfdMakeFigure14()
%
% AUTHORS. TITLE. JOURNAL. YEAR.
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
meshData([1 2 4])     = dfdMakeFigure14AcrossSubjects(whichSubjectsRaw,figureDir,dataDir,saveFigures,threshold);


%% Get denoising results for tsss data
whichSubjectsTSSS   = [14,16,18,20];       % TSSS
meshData([1 3 5])     = dfdMakeFigure14AcrossSubjects(whichSubjectsTSSS,figureDir,dataDir,saveFigures,threshold);

%% Make bar graph

% Get top ten across conditions for each subject
snr_diff = zeros(length(whichSubjectsRaw),size(meshData,2)-1,3); % All broadband columns
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
            snr_pre   = meshData{2}(icond,finalpcchan(whichSubject,:),whichSubject); % Second column is undenoised Broadband
            snr_post  = meshData{thisColumn}(icond,finalpcchan(whichSubject,:),whichSubject);
            
            snr_diff(whichSubject,thisColumn-1,icond) = mean(snr_post-snr_pre);
        end
    end
end

%% Plot figure
fH = figure('position',[0,300,700,300]);
% define what the different conditions are
types = {'MEG Raw - Sanity check','TSSS','MEG Denoise','MEG Denoise + TSSS'}; %
% re-arrange the order of the bars
neworder = [1:4];
newtypes = types(neworder);
colors = dfdGetColors(3);

snr_diff2 = snr_diff(:,neworder,:);
nnull = length(neworder);
for icond = 1:3
    subplot(1,3,icond);
    % mean and sem across subjects
    mn  = mean(snr_diff2(:,:,icond));
    sem = std(snr_diff2(:,:,icond))/sqrt(4);
    bar(1:nnull, mn,'EdgeColor','none','facecolor',colors(icond,:)); hold on
    errorbar2(1:nnull,mn,sem,1,'-','color',colors(icond,:));
    % format figure and make things pretty
    set(gca,'xlim',[0.2,nnull+0.8],'ylim',[-1,5]);
    makeprettyaxes(gca,9,9);
    set(gca,'XTickLabel',types(neworder));
    set(gca,'XTickLabelRotation',45);
%     set(get(gca,'XLabel'),'Rotation',45);
    ylabel('Difference in broadband SNR (post-pre)')
end

if saveFigures
    figurewrite(fullfile(figureDir,'figure14b_bargraphTSSS'),[],0,'.',1);
end