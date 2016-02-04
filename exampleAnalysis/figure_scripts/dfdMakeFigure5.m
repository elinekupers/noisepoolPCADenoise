function dfdMakeFigure5()
%% Function to reproduce Figure 5 (Spatialmap) from example subject 
%
% dfdMakeFigure5()
%
% AUTHORS. TITLE. JOURNAL. YEAR.
%
% This figure will show an interpolated spatial map of the SNR values in 
% each channel for the stimulus locked signal, broadband signals before and
% after using the denoising algorithm. The three separate conditions (Full,
% left, right hemifield stimulation are shown separately). 
%
% This function assumes that data is downloaded with the DFDdownloaddata
% function. 

%% Choices to make:                                              
whichSubject    = 9;        % Subject 1 is the example subject.
figureDir       = fullfile(dfdRootPath, 'exampleAnalysis', 'figures_rm1epoch'); % Where to save images?
dataDir         = fullfile(dfdRootPath, 'exampleAnalysis', 'data');    % Where to save data?
saveFigures     = true;     % Save figures in the figure folder?
threshold       = 2;        % Set threshold for colormap. If no threshold set value to 0

if whichSubject > 8;
    data_hdr        = 'neuromag360_sample_hdr_combined.mat';
    cfg             = 'neuromag360_sample_cfg_combined.mat';
else
    data_hdr        = [];
end

%% Load denoised data of example subject
[data] = prepareData(dataDir,whichSubject,5);
bb = data{1};
sl = data{2};

%% Plot stimulus-locked signal, broadband before and after denoising on sensormap
figure('position',[1,600,1400,800]);
condNames = {'Stim Full','Stim Left','Stim Right'};
for icond = 1:3
    % get stimulus-locked snr
    sl_snr1 = getsignalnoise(sl.results.origmodel(1),icond, 'SNR');
    %clims_sl = [0, max(sl_snr1)];
    clims_sl = [0,25.6723];
    % get broadband snr for before and after denoising
    ab_snr1 = getsignalnoise(bb.results.origmodel(1),  icond, 'SNR');
    ab_snr2 = getsignalnoise(bb.results.finalmodel(1), icond, 'SNR');
    clims_ab = [0, max([ab_snr1, 8.4445])];
    clims_ab = [0, 7.4445];
    %clims_ab = [0, max([ab_snr1, ab_snr2])];
    
    if whichSubject < 9; % Assuming NeuroMag360 data has no badChannels (CHECK!)
        % convert back into 157-channel space
        ab_snr1 = to157chan(ab_snr1,~bb.badChannels,'nans');
        ab_snr2 = to157chan(ab_snr2,~bb.badChannels,'nans');
        sl_snr1 = to157chan(sl_snr1,~sl.badChannels,'nans');
    end
    
    % Set threshold
    ab_snr1(abs(ab_snr1) < threshold) = 0;
    ab_snr2(abs(ab_snr2) < threshold) = 0;
    sl_snr1(abs(sl_snr1) < threshold) = 0;
 
    % plot spatial maps
    subplot(3,3,(icond-1)*3+1)
    [~,ch] = megPlotMap(sl_snr1,clims_sl,gcf,'parula',sprintf('%s : Stimulus Locked Original', condNames{icond}),data_hdr,cfg);
    makeprettyaxes(gca,9,9);
    makeprettyaxes(ch,9,9);
    title(sprintf('SL no DN %s', condNames{icond}))
    
    subplot(3,3,(icond-1)*3+2)
    [~,ch] = megPlotMap(ab_snr1,clims_ab,gcf,'parula',sprintf('%s Original', condNames{icond}),data_hdr,cfg);
    makeprettyaxes(gca,9,9);
    makeprettyaxes(ch,9,9);
    title(sprintf('Broadband Pre %s', condNames{icond}))
    
    subplot(3,3,(icond-1)*3+3)
    [~,ch] = megPlotMap(ab_snr2,clims_ab,gcf,'parula',sprintf('%s : Denoised PC %d',condNames{icond}, bb.results.pcnum(1)),data_hdr,cfg);
    makeprettyaxes(gca,9,9);
    makeprettyaxes(ch,9,9);
    title(sprintf('Broadband Post %s', condNames{icond}))
end

if saveFigures
    figurewrite(sprintf(fullfile(figureDir,'figure5_examplesubject%d_parula_thresh%d_varThresh3'),whichSubject, threshold),[],0,'.',1);
end

