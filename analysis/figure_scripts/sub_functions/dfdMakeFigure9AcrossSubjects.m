function dfdMakeFigure9AcrossSubjects()

%% Function to reproduce Figure 9 (Spatialmap) across all subjects
%
% dfdMakeFigure9AcrossSubjects()
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
whichSubjects    = [1:8]; %[13,15,17,19];        % Subject 1 is the example subject.
% whichSubjects    = [9:12];%[14,16,18,20];        % Subject 1 is the example subject.
figureDir       = fullfile(dfdRootPath, 'analysis', 'figures'); % Where to save images?
dataDir         = fullfile(dfdRootPath, 'analysis', 'data');    % Where to save data?
saveFigures     = true;     % Save figures in the figure folder?
threshold       = 0;

%% Compute SNR across subjects
contrasts = [1 0 0; 0 1 0; 0 0 1; 0 1 -1]; % Full, Left, Right and L-R

computeSNR    = @(x) nanmean(x,3) ./ nanstd(x, [], 3);
% computeSignal = @(x) nanmean(x,3);


contrastNames = {
    'Full'...
    'Left'...
    'Right'...
    'Left-Right'
    };

%% Load denoised data of all subjects

for whichSubject = whichSubjects
    subjnum = find(whichSubjects==whichSubject);
    data = prepareData(dataDir,whichSubject,9);
    bb(subjnum) = data;
    
    num_channels = size(bb(subjnum).results.origmodel.beta,2);
    num_boots    = size(bb(subjnum).results.origmodel.beta,3);
    num_contrasts = length(contrasts);
    
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
    
    sBBBeforeAcrossSubjects(:,:,subjnum) = to157chan(sBBBefore', ~bb(subjnum).badChannels,'nans');
    sBBAfterAcrossSubjects(:,:,subjnum) = to157chan(sBBAfter', ~bb(subjnum).badChannels,'nans');
        
end



%% Plot stimulus-locked signal, broadband before and after denoising on sensormap
figure('position',[1,600,1400,800], 'Name', 'Figure 9, group data', 'NumberTitle', 'off');
for icond = 1:numel(contrastNames)
    
    % get broadband snr before denoising
    %     ab_snr1 = nanmean(sBBBeforeAcrossSubjects,3)  ./ sem(sBBBeforeAcrossSubjects,3);
    ab_snr1 = nanmean(sBBBeforeAcrossSubjects,3);
    
    % threshold
    ab_snr1(abs(ab_snr1) < threshold) = 0;
    
    % get broadband snr after denoising
    %     ab_snr2 = nanmean(sBBAfterAcrossSubjects,3)  ./ sem(sBBAfterAcrossSubjects,3);
    ab_snr2 = nanmean(sBBAfterAcrossSubjects,3);
    
    % threshold
    ab_snr2(abs(ab_snr2) < threshold) = 0;
    
    
    % Define ranges colormap
%     clims_ab = [-6.4445,6.4445];
    if max(unique(whichSubjects)) < 9
        clims_ab = [-8,8];
        if icond == 4; clims_ab = 3.5363 * [-1 1];  end;
    else
        clims_ab = [-4,4];
        if icond == 4; clims_ab = [-4, 4];  end;
    end
    cmap = 'bipolar';
    
    
    if size(ab_snr1,2) > 157
        % Combine channels
        ab_snr1 = dfd204to102(ab_snr1(icond,:));
        ab_snr2 = dfd204to102(ab_snr2(icond,:));
    else
        ab_snr1 = ab_snr1(icond,:);
        ab_snr2 = ab_snr2(icond,:);
    end
    
    % plot spatial maps   
    subplot(4,2,(icond-1)*2+1)
    [~,ch] = megPlotMap((ab_snr1),clims_ab,gcf,cmap,...
        sprintf('%s Original', contrastNames{icond}));
    makeprettyaxes(ch,9,9); 
    if icond == 4; ticks = [-5,-2.5,0,2.5,5]; else ticks = [-8,-4,0,4,8]; end; set(ch,'YTick',ticks);
    title(sprintf('Broadband Pre %s', contrastNames{icond}))
    
    subplot(4,2,(icond-1)*2+2)
    [~,ch] = megPlotMap((ab_snr2),clims_ab,gcf,cmap,...
        sprintf('%s : Denoised PC %d',contrastNames{icond}, bb(end).results.pcnum(1)));
    makeprettyaxes(ch,9,9);
    if icond == 4; ticks = [-5,-2.5,0,2.5,5]; else ticks = [-8,-4,0,4,8]; end; set(ch,'YTick',ticks);
    title(sprintf('Broadband Post %s', contrastNames{icond}))
    drawnow();
   
end


if saveFigures
     hgexport(gcf,fullfile(figureDir, sprintf('figure9_AcrossSubject%d_bipolar_threshold%d',whichSubject, threshold)));
%    figurewrite(fullfile(figureDir, sprintf('figure9_AcrossSubject%d_bipolar_threshold%d',whichSubject, threshold)),[],0,'.',1);
end
