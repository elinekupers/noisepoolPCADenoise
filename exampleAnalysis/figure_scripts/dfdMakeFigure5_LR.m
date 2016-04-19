function dfdMakeFigure5_LR()

%% Choices to make:                                              
whichSubject    = 8;        % Subject 1 is the example subject. NB: Can only run one subject at the time.
figureDir       = fullfile(dfdRootPath, 'exampleAnalysis', 'figures_rm1epoch'); % Where to save images?
dataDir         = fullfile(dfdRootPath, 'exampleAnalysis', 'data');    % Where to save data?
saveFigures     = true;     % Save figures in the figure folder?
threshold       = 0;

%% Load denoised data of example subject
[data] = prepareData(dataDir,whichSubject,5);
bb = data{1};
sl = data{2};

%% Plot stimulus-locked signal, broadband before and after denoising on sensormap
figure('position',[1,600,1400,800]);
condNames = {'Stim Full','Stim Left','Stim Right'};

% get stimulus-locked snr
sl_snr1_L = getsignalnoise(sl.results.origmodel(1), 2, 'SNR');
sl_snr1_R = getsignalnoise(sl.results.origmodel(1), 3, 'SNR');


%clims_sl = [0, max(sl_snr1)];
clims_sl = [-25.6723,25.6723];

ab_snr1_L = getsignalnoise(bb.results.origmodel(1),  2, 'SNR');
ab_snr1_R = getsignalnoise(bb.results.origmodel(1),  3, 'SNR');

ab_snr2_L = getsignalnoise(bb.results.finalmodel(1), 2, 'SNR');
ab_snr2_R = getsignalnoise(bb.results.finalmodel(1), 3, 'SNR');


clims_ab = [- 8.4445, 8.4445];  
%clims_ab = [0, max([ab_snr1, ab_snr2])];

% convert back into 157-channel space
ab_snr1a_LmnR = to157chan(ab_snr1_L,~bb.badChannels,'nans') - to157chan(ab_snr1_R,~bb.badChannels,'nans');
ab_snr2a_LmnR = to157chan(ab_snr2_L,~bb.badChannels,'nans') - to157chan(ab_snr2_R,~bb.badChannels,'nans');
sl_snr1a_LmnR = to157chan(sl_snr1_L,~bb.badChannels,'nans') - to157chan(sl_snr1_R,~bb.badChannels,'nans');
   

% Threshold
ab_snr1a_LmnR(abs(ab_snr1a_LmnR) < threshold) = 0;
ab_snr2a_LmnR(abs(ab_snr2a_LmnR) < threshold) = 0;
sl_snr1a_LmnR(abs(sl_snr1a_LmnR) < threshold) = 0;

% plot spatial maps
subplot(311)
[~,ch] = megPlotMap(sl_snr1a_LmnR,clims_sl,gcf,'bipolar',sprintf('%s : Stimulus Locked Original', 'Left minus Right'));
makeprettyaxes(gca,9,9);
makeprettyaxes(ch,9,9);
title(sprintf('SL no DN %s', 'Left minus Right'))

subplot(312)
[~,ch] = megPlotMap(ab_snr1a_LmnR,clims_ab,gcf,'bipolar',sprintf('%s Original', 'Left minus Right'));
makeprettyaxes(gca,9,9);
makeprettyaxes(ch,9,9);
title(sprintf('Broadband Pre %s', 'Left minus Right'))

subplot(313)
[~,ch] = megPlotMap(ab_snr2a_LmnR,clims_ab,gcf,'bipolar',sprintf('%s : Denoised PC %d', 'Left minus Right', bb.results.pcnum(1)));
makeprettyaxes(gca,9,9);
makeprettyaxes(ch,9,9);
title(sprintf('Broadband Post %s', 'Left minus Right'))

if saveFigures
    printnice(gcf, 0, figureDir, sprintf('figure5_LR_examplesubject%d_bipolar_threshold%d',whichSubject, threshold));
end

end