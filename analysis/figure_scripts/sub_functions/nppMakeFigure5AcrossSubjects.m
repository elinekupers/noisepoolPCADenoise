function nppMakeFigure5AcrossSubjects

%% Function to reproduce Figure 5 (Spatialmap) across all subjects
%
% nppMakeFigure5AcrossSubjects()
%
% AUTHORS. TITLE. JOURNAL. YEAR.
%
% This figure will show an interpolated spatial map of the SNR values in
% each channel for the stimulus locked signal, broadband signals before
% using the denoising algorithm. The three separate conditions (Full,
% left, right hemifield stimulation are shown separately).
%
% This function assumes that data is downloaded with the nppDownloaddata
% function.

%% Choices to make:
whichSubjects    = 1:8;        % Subject 1 is the example subject.
% whichSubjects    = [9:12];%[14,16,18,20];        % Subject 1 is the example subject.
figureDir       = fullfile(nppRootPath, 'analysis', 'figures'); % Where to save images?
dataDir         = fullfile(nppRootPath, 'analysis', 'data');    % Where to save data?
saveFigures     = true;     % Save figures in the figure folder?
threshold       = 0;
climsSL         = [-25.6723,25.6723];
climsAB         = [-8,8];
cmap            = 'bipolar';
yscaleAB        = [repmat([-8,-4,0,4,8],3,1);[-5,-2.5,0,2.5,5]];

%% Compute SNR across subjects
contrasts = [eye(3); 0 1 -1];
contrasts = bsxfun(@rdivide, contrasts, sqrt(sum(contrasts.^2,2))); % Full, Left, Right and L-R
computeSNR    = @(x) nanmean(x,3) ./ nanstd(x, [], 3);

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
    bb = data{1}; sl = data{2}; clear data;
    
    % Get information, assuming this is the same for SL and BB data
    numChannels  = size(sl.results.origmodel.beta,2);
    numBoots     = size(sl.results.origmodel.beta,3);
    numContrasts = length(contrasts);
    
    % Stimulus-locked: Compute SNR for contrasts
    tmp_data = reshape(sl.results.origmodel.beta,size(contrasts,2),[]);
    tmp      = contrasts*tmp_data;
    tmp      = reshape(tmp, numContrasts, numChannels, numBoots);
    sSL      = computeSNR(tmp)';
    
    % Broadband before denoising: Compute SNR for contrasts
    tmp_data = reshape(bb.results.origmodel.beta,size(contrasts,2),[]);
    tmp = contrasts*tmp_data;
    tmp = reshape(tmp, numContrasts, numChannels,numBoots);
    sBBBefore = computeSNR(tmp)';
    
    % Prepare array
    if subjnum == 1
        sSLAcrossSubjects = NaN(size(contrasts,1),length(sl(1).badChannels), length(whichSubjects));
        sBBBeforeAcrossSubjects = sSLAcrossSubjects;
    end
    
    % Update array with data converted to channel space
    sSLAcrossSubjects(:,:,subjnum) = to157chan(sSL', ~sl.badChannels,'nans');
    sBBBeforeAcrossSubjects(:,:,subjnum) = to157chan(sBBBefore', ~bb.badChannels,'nans');
    
end

%% Plot stimulus-locked signal, broadband before on sensormap
figure('position',[1,600,1400,800]); set(gcf, 'Name', 'Figure 5, Group Data', 'NumberTitle', 'off')

for icond = 1:numel(contrastNames)
    
    % get stimulus-locked snr
    sl_snr1 = nanmean(sSLAcrossSubjects,3);
    
    % threshold
    sl_snr1(abs(sl_snr1) < threshold) = 0;
    
    % get broadband snr before denoising
    ab_snr1 = nanmean(sBBBeforeAcrossSubjects,3);
    
    % threshold
    ab_snr1(abs(ab_snr1) < threshold) = 0;
    
    % Change ranges colormap when contrast is Left-Right
    if icond == 4; climsAB = [-3.5363,3.5363];  end;
   
    % Make sure the channels are combined if CiNet data
    if size(sl_snr1,2) > 157
        % Combine channels
        sl_snr1 = npp204to102(sl_snr1(icond,:));
        ab_snr1 = npp204to102(ab_snr1(icond,:));
    else
        sl_snr1 = sl_snr1(icond,:);
        ab_snr1 = ab_snr1(icond,:);
    end
    
    % Plot spatial maps
    subplot(4,2,(icond-1)*2+1)
    [~,ch] = megPlotMap((sl_snr1),climsSL,gcf,cmap, sprintf('%s : Stimulus-locked Original', contrastNames{icond}));
    makeprettyaxes(ch,9,9);   
    
    subplot(4,2,(icond-1)*2+2)
    [~,ch] = megPlotMap((ab_snr1),climsAB,gcf,cmap, sprintf('%s : Broadband Original', contrastNames{icond}));
    makeprettyaxes(ch,9,9);
    set(ch,'YTick',yscaleAB(icond,:));

    
end

if saveFigures
%      figurewrite(fullfile(figureDir, sprintf('figure5_AcrossSubjects%d_threshold%d',whichSubject, threshold)),[],0,'.',1);
     hgexport(gcf,fullfile(figureDir, sprintf('figure5_AcrossSubjects%d_threshold%d',whichSubject, threshold)));
end
