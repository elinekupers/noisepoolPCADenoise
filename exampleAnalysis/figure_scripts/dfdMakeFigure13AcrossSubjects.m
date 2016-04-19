function dfdMakeFigure13AcrossSubjects
%
% NEEDS SOME CLEANING


%% Function to reproduce Figure 5 (Spatialmap) across all subjects
%
% dfdMakeFigure5AcrossSubjects()
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
whichSubjects    = [14,16,18,20];        % Subject 1 is the example subject.
% whichSubjects    = 9:12;        % Subject 1 is the example subject.
figureDir       = fullfile(dfdRootPath, 'exampleAnalysis', 'figures_rm1epoch_CiNet'); % Where to save images?
dataDir         = fullfile(dfdRootPath, 'exampleAnalysis', 'data');    % Where to save data?
saveFigures     = true;     % Save figures in the figure folder?
threshold       = 2;

%% Compute SNR across subjects
contrasts = [1 0 0; 0 1 0; 0 0 1; 0 1 -1]; % Full, Left, Right and L-R

computeSNR    = @(x) nanmean(x,3) ./ nanstd(x, [], 3);
computeSignal = @(x) nanmean(x,3);


contrastNames = {
    'Full'...
    'Left'...
    'Right'...
    'Left-Right'
    };

%% Load denoised data of all subjects
 
for whichSubject = whichSubjects
    subjnum = find(whichSubjects==whichSubject);
    data = prepareData(dataDir,whichSubject,5);
    bb(subjnum) = data{1};
    sl(subjnum) = data{2};
    
    % SL
    num_channels = size(sl(subjnum).results.origmodel.beta,2);
    num_boots    = size(sl(subjnum).results.origmodel.beta,3);
    num_contrasts = length(contrasts);
    
    tmp_data = reshape(sl(subjnum).results.origmodel.beta,3,[]);
    tmp = contrasts*tmp_data;
    tmp = reshape(tmp, num_contrasts, num_channels, num_boots);
    sSL = computeSNR(tmp)';
    
    % BB before
    tmp_data = reshape(bb(subjnum).results.origmodel.beta,3,[]);
    tmp = contrasts*tmp_data;
    tmp = reshape(tmp, num_contrasts, num_channels,num_boots);
    sBBBefore = computeSNR(tmp)';
    
    % BB before
    tmp_data = reshape(bb(subjnum).results.finalmodel.beta,3,[]);
    tmp = contrasts*tmp_data;
    tmp = reshape(tmp, num_contrasts, num_channels,num_boots);
    sBBAfter = computeSNR(tmp)';
    
    if subjnum == 1
        sSLAcrossSubjects = NaN(size(contrasts,1),length(sl(1).badChannels), length(whichSubjects));
        sBBBeforeAcrossSubjects = sSLAcrossSubjects;
        sBBAfterAcrossSubjects  = sSLAcrossSubjects;
    end
    
    sSLAcrossSubjects(:,:,subjnum) = to157chan(sSL', ~sl(subjnum).badChannels,'nans');
    sBBBeforeAcrossSubjects(:,:,subjnum) = to157chan(sBBBefore', ~bb(subjnum).badChannels,'nans');
    sBBAfterAcrossSubjects(:,:,subjnum) = to157chan(sBBAfter', ~bb(subjnum).badChannels,'nans');

    
end



%% Plot stimulus-locked signal, broadband before and after denoising on sensormap
figure('position',[1,600,1400,800]);
condNames = {'Stim Full','Stim Left','Stim Right'};
n = 0;



% stimulus locked
data{1} = sSLAcrossSubjects;
data{2} = sBBBeforeAcrossSubjects;
data{3} = sBBAfterAcrossSubjects;
cmap = bipolar;
str = {'SL-pre' 'BB-pre' 'BB-post'};

figure,set(gcf, 'Name', 'RAW')
for row = 1:4 % stimulus contrasts
    for col = 1:3 % types of analyses (sl, bb-pre, bb-post)
        subplot(4,3,3*(row-1)+col),
        if col == 1, clim = [-15 15]; else clim = [-4 4]; end
        megPlotMap(dfd204to102(squeeze(mean(data{col}(row,:,:),3))), ...
            clim, [], cmap);        
        if row == 1, title(str{col}); end
    end
end

if saveFigures
    figurewrite(sprintf(fullfile(figureDir,'figure13_AcrossSubject%d_bipolar_threshold%d_tsss'),whichSubject, threshold),[],0,'.',1);
end
