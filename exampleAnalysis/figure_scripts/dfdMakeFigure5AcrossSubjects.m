function dfdMakeFigure5AcrossSubjects()

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
whichSubjects    = 1:8;        % Subject 1 is the example subject.
figureDir       = fullfile(dfdRootPath, 'exampleAnalysis', 'figures_rm1epoch'); % Where to save images?
dataDir         = fullfile(dfdRootPath, 'exampleAnalysis', 'data');    % Where to save data?
saveFigures     = true;     % Save figures in the figure folder?

%% Compute SNR across subjects
contrasts = [1 0 0; 0 1 0; 0 0 1; 0 1 -1]; % Full, Left, Right and L-R
computeSNR = @(x) nanmean(x,3) ./ nanstd(x, [], 3);

contrastNames = {
    'Full'...
    'Left'...
    'Right'...
    'Left-Right'
    };

%% Load denoised data of all subjects

for whichSubject = whichSubjects
    
    data = prepareData(dataDir,whichSubject,5);
    bb(whichSubject,:) = data{1};
    sl(whichSubject,:) = data{2};
    
    % SL
    num_channels = size(sl(whichSubject).results.origmodel.beta,2);
    num_boots    = size(sl(whichSubject).results.origmodel.beta,3);
    num_contrasts = length(contrasts);
    
    tmp_data = reshape(sl(whichSubject).results.origmodel.beta,3,[]);
    tmp = contrasts*tmp_data;
    tmp = reshape(tmp, num_contrasts, num_channels, num_boots);
    snrSL = computeSNR(tmp)';
    
    % BB before
    tmp_data = reshape(bb(whichSubject).results.origmodel.beta,3,[]);
    tmp = contrasts*tmp_data;
    tmp = reshape(tmp, num_contrasts, num_channels,num_boots);
    snrBB = computeSNR(tmp)';
    
    
    snrSLAcrossSubjects(:,:,whichSubject) = to157chan(snrSL', ~sl(whichSubject).badChannels,'nans');
    snrBBAcrossSubjects(:,:,whichSubject) = to157chan(snrBB', ~bb(whichSubject).badChannels,'nans');

end



%% Plot stimulus-locked signal, broadband before and after denoising on sensormap
figure('position',[1,600,1400,800]);
condNames = {'Stim Full','Stim Left','Stim Right'};
for icond = 1:length(contrasts)
    
    % get stimulus-locked snr
    sl_snr1 = nanmean(snrSLAcrossSubjects,3);
    
    % get broadband snr for before
    ab_snr1 = nanmean(snrBBAcrossSubjects,3);
    %clims_sl = [0, max(sl_snr1)];
    clims_sl = [0,25.6723];
    clims_sl = [0,20];
%     % get broadband snr for before and after denoising
%     ab_snr1 = getsignalnoise(bb.results.origmodel(1),  icond, 'SNR');
%     ab_snr2 = getsignalnoise(bb.results.finalmodel(1), icond, 'SNR');
%     clims_ab = [0, 10.4445];
	clims_ab = [0,8];
    %clims_ab = [0, max([ab_snr1, ab_snr2])];
    
%     % convert back into 157-channel space
%     ab_snr1a = to157chan(ab_snr1,~bb.badChannels,'nans');
%     ab_snr2a = to157chan(ab_snr2,~bb.badChannels,'nans');
%     sl_snr1a = to157chan(sl_snr1,~sl.badChannels,'nans');
%     
    % plot spatial maps
    subplot(4,2,(icond-1)*2+1)
    [~,ch] = megPlotMap(sl_snr1(icond,:),clims_sl,gcf,'jet',sprintf('%s : Stimulus Locked Original', contrastNames{icond}));
    makeprettyaxes(gca,9,9);
    makeprettyaxes(ch,9,9);
    title(sprintf('SL no DN %s', contrastNames{icond}))
    
    subplot(4,2,(icond-1)*2+2)
    [~,ch] = megPlotMap(ab_snr1(icond,:),clims_ab,gcf,'jet',sprintf('%s Original', contrastNames{icond}));
    makeprettyaxes(gca,9,9);
    makeprettyaxes(ch,9,9);
    title(sprintf('Broadband Pre %s', contrastNames{icond}))
    
%     subplot(3,3,(icond-1)*3+3)
%     [~,ch] = megPlotMap(ab_snr2a,clims_ab,gcf,'jet',sprintf('%s : Denoised PC %d',condNames{icond}, bb.results.pcnum(1)));
%     makeprettyaxes(gca,9,9);
%     makeprettyaxes(ch,9,9);
%     title(sprintf('Broadband Post %s', condNames{icond}))
end

if saveFigures
    figurewrite(sprintf(fullfile(figureDir,'figure5_AcrossSubject%d'),whichSubject),[],0,'.',1);
end
