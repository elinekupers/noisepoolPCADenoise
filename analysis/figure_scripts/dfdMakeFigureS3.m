function dfdMakeFigure3()
%% Function to reproduce Supplementary Figure 3 
%
% dfdMakeFigureS3()
%
% Eline Kupers, Helena X. Wang, Kaoru Amano, Kendrick N. Kay, David J.
% Heeger, Jonathan Winawer. (YEAR) Broadband spectral responses in visual
% cortex revealed by a new MEG denoising algorithm.
% (JOURNAL. VOLUME. ISSUE. DOI.)
%
% This figure will show an interpolated spatial map of the SNR values in
% each channel for the stimulus locked signal, broadband signals before and
% after using the denoising algorithm for all individual subjects.
% The three separate conditions (Full, left, right hemifield stimulation are
% shown separately).
%
% This function assumes that data is downloaded with the dfdDownloadsampledata
% function.

%% Choices to make:

figureDir       = fullfile(dfdRootPath, 'analysis', 'figures'); % Where to save images?
dataDir         = fullfile(dfdRootPath, 'analysis', 'data');    % Where to save data?
saveFigures     = true;     % Save figures in the figure folder?
threshold       = 0;        % Set threshold for colormap. If no threshold set value to 0
cfg             = [];
data_hdr        = [];

%% Plot stimulus-locked signal, broadband after denoising on sensormap
fH1 = figure('position',[1,600,1400,800]); set(gcf, 'Name', 'SF2, Individual subjects BB Post denoising', 'NumberTitle', 'off');

for whichSubject = 1:8
    %% Load denoised data of example subject
    bb = prepareData(dataDir,whichSubject,8);
    
    %% Plot stimulus-locked signal, broadband after denoising on sensormap
    contrasts = [eye(3); [0 1 -1]/sqrt(2)];
    for icond = 1:4
        
        % get broadband snr for after denoising
        ab_snr2 = getsignalnoise(bb.results.finalmodel(1), contrasts(icond,:), 'SNR',bb.badChannels);
        
        
        % Replace bad channels with NaNs
        ab_snr2 = to157chan(ab_snr2,~bb.badChannels,'nans');
        
        if length(bb.badChannels)>157
            ab_snr2 = dfd204to102(ab_snr2);
        end
        
        % Threshold
        ab_snr2(abs(ab_snr2) < threshold) = 0;
        
        % Set colormap limits
        if icond <= 3, clims_ab = 8.4445 * [-1 1]; else, clims_ab = 3.5363 * [-1,1]; end
        
        subplot(4,8,whichSubject+(icond*8)-8)
        [~,ch] = megPlotMap(ab_snr2,clims_ab,gcf,'bipolar',sprintf('S%d',whichSubject),data_hdr,cfg);
        makeprettyaxes(ch,9,9);
        set(ch,'YTick',[-8,-4,0,4,8]); if icond > 3; set(ch,'YTick',[-5,-2.5,0,2.5,5]); end
        drawnow();
        
    end
    
end

if saveFigures
    % Figurewrite takes a long time saving, so let's use hgexport for now
%     figurewrite(fullfile(figureDir, sprintf('SF2_individualsubject_thresh%d', threshold)),[],0,'.',1);
    hgexport(fH1, fullfile(figureDir, sprintf('S3_individualsubject_thresh%d', threshold)));
end

