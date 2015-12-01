function dfdMakeFigure8AcrossSubjects()

%% Function to reproduce Figure 8 (Spatialmap) across all subjects
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
    
    % SL
    num_channels = size(bb(whichSubject).results.origmodel.beta,2);
    num_boots    = size(bb(whichSubject).results.origmodel.beta,3);
    num_contrasts = length(contrasts);
    
    tmp_data = reshape(bb(whichSubject).results.origmodel.beta,3,[]);
    tmp = contrasts*tmp_data;
    tmp = reshape(tmp, num_contrasts, num_channels, num_boots);
    snrBBB = computeSNR(tmp)';
    
    % BB before
    tmp_data = reshape(bb(whichSubject).results.finalmodel.beta,3,[]);
    tmp = contrasts*tmp_data;
    tmp = reshape(tmp, num_contrasts, num_channels,num_boots);
    snrBBA = computeSNR(tmp)';
    
    
    snrBBBAcrossSubjects(:,:,whichSubject) = to157chan(snrBBB', ~bb(whichSubject).badChannels,'nans');
    snrBBAAcrossSubjects(:,:,whichSubject) = to157chan(snrBBA', ~bb(whichSubject).badChannels,'nans');
    
end



%% Plot stimulus-locked signal, broadband before and after denoising on sensormap
figure('position',[1,600,1400,800]);
condNames = {'Stim Full','Stim Left','Stim Right'};
for icond = 1:length(contrasts)
    
    % get stimulus-locked snr
    ab_snr1 = nanmean(snrBBBAcrossSubjects,3);
    
    % get broadband snr for before
    ab_snr2 = nanmean(snrBBAAcrossSubjects,3);
    
    if icond == 4;
        clims_sl = [-10,10];
        clims_ab = [-4,4];
    else
        clims_sl = [0,25.6723];
        clims_sl = [0,20];
        
        %     clims_ab = [0, 10.4445];
        clims_ab = [0,8];
    end
    
    % plot spatial maps
%     subplot(4,2,(icond-1)*2+1)
%     [~,ch] = megPlotMap(sl_snr1(icond,:),clims_sl,gcf,'jet',sprintf('%s : Stimulus Locked Original', contrastNames{icond}));
%     makeprettyaxes(gca,9,9);
%     makeprettyaxes(ch,9,9);
%     title(sprintf('SL no DN %s', contrastNames{icond}))
    
    subplot(4,2,(icond-1)*2+1)
    [~,ch] = megPlotMap(ab_snr1(icond,:),clims_ab,gcf,'jet',sprintf('%s Original', contrastNames{icond}));
    makeprettyaxes(gca,9,9);
    makeprettyaxes(ch,9,9);
    title(sprintf('Broadband Pre %s', contrastNames{icond}))
    
    subplot(4,2,(icond-1)*2+2)
    [~,ch] = megPlotMap(ab_snr2(icond,:),clims_ab,gcf,'jet',sprintf('%s : Denoised PC %02d',contrastNames{icond}, bb(1).results.pcnum(1)));
    makeprettyaxes(gca,9,9);
    makeprettyaxes(ch,9,9);
    title(sprintf('Broadband Post %s', contrastNames{icond}))
end

if saveFigures
    figurewrite(sprintf(fullfile(figureDir,'figure8_AcrossSubject%d'),whichSubject),[],0,'.',1);
end
