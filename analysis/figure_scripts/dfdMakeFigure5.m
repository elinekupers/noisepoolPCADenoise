function dfdMakeFigure5(whichSubject)
%% Function to reproduce Figure 5 (Spatialmap) from example subject 
%
% dfdMakeFigure5()
%
% Eline Kupers, Helena X. Wang, Kaoru Amano, Kendrick N. Kay, David J.
% Heeger, Jonathan Winawer. (YEAR) Broadband spectral responses in visual
% cortex revealed by a new MEG denoising algorithm.
% (JOURNAL. VOLUME. ISSUE. DOI.)
%
% This figure will show an interpolated spatial map of the SNR values in 
% each channel for the stimulus locked signal, broadband signals before and
% after using the denoising algorithm. The three separate conditions (Full,
% left, right hemifield stimulation are shown separately). 
%
% This function assumes that data is downloaded with the dfdDownloadsampledata
% function. 

%% Choices to make:                                              
if notDefined('whichSubject')
    whichSubject    = 1;        % Subject 1 is the example subject.
end
figureDir       = fullfile(dfdRootPath, 'analysis', 'figures'); % Where to save images?
dataDir         = fullfile(dfdRootPath, 'analysis', 'data');    % Where to save data?
saveFigures     = true;     % Save figures in the figure folder?
threshold       = 0;        % Set threshold for colormap. If no threshold set value to 0
cfg             = [];
data_hdr        = [];


%% Load denoised data of example subject
[data] = prepareData(dataDir,whichSubject,5);
bb = data{1};
sl = data{2};

%% Plot stimulus-locked signal, broadband before and after denoising on sensormap
figure('position',[1,600,1400,800]); set(gcf, 'Name', 'Figure 5, Example subject', 'NumberTitle', 'off');
contrastNames = {'Stim Full','Stim Left','Stim Right','Left minus Right'};
contrasts = [eye(3); 0 1 -1];
contrasts = bsxfun(@rdivide, contrasts, sqrt(sum(contrasts.^2,2)));
yscaleAB = [repmat([-8,-4,0,4,8],3,1);[-5,-2.5,0,2.5,5]];
climsSL = [-25.6723,25.6723];
climsAB = [-8.4445, 8.4445];

for icond = 1:size(contrasts,1)
    % get stimulus-locked snr
    sl_snr1 = getsignalnoise(sl.results.origmodel(1),contrasts(icond,:), 'SNR',sl.badChannels);
    % get broadband snr for before and after denoising
    ab_snr1 = getsignalnoise(bb.results.origmodel(1),  contrasts(icond,:), 'SNR',bb.badChannels);
    


    sl_snr1 = to157chan(sl_snr1,~sl.badChannels,'nans');
    ab_snr1 = to157chan(ab_snr1,~bb.badChannels,'nans');

    if length(bb.badChannels)>157 && length(bb.badChannels) <= 204
        sl_snr1 = dfd204to102(sl_snr1);
        ab_snr1 = dfd204to102(ab_snr1);
    end  
    
    % Threshold if requested
    ab_snr1(abs(ab_snr1) < threshold) = 0;
    sl_snr1(abs(sl_snr1) < threshold) = 0;
   
     if icond > 3 % then we are plotting l-r rather than one condition and change the colormap limits
        climsAB = [-3.5363, 3.5363]; 
     end

    % plot spatial maps
    subplot(4,2,(icond-1)*2+1)
    [~,ch] = megPlotMap(sl_snr1,climsSL,gcf,'bipolar',sprintf('%s : Stimulus Locked Original', contrastNames{icond}),data_hdr,cfg);
    makeprettyaxes(ch,9,9);
    set(ch,'YTick',[-20,-10,0,10,20]);
    
    subplot(4,2,(icond-1)*2+2)
    [~,ch] = megPlotMap(ab_snr1,climsAB,gcf,'bipolar',sprintf('%s : Broadband Original', contrastNames{icond}),data_hdr,cfg);
    makeprettyaxes(ch,9,9);
    set(ch,'YTick',yscaleAB(icond,:));

end

 if saveFigures
    % Note: our function figurewrite is extremely slow with Matlab 2016b,
    % therefore we use hgexport()
%      figurewrite(fullfile(figureDir, sprintf('figure5_examplesubject%d_thresh%d',whichSubject, threshold)),[],0,'.',1);
    hgexport(gcf,fullfile(figureDir, sprintf('figure5_examplesubject%d_thresh%d',whichSubject, threshold)))
 end

%% Now call dfdMakeFigure3AcrossSubjects
dfdMakeFigure5AcrossSubjects();


