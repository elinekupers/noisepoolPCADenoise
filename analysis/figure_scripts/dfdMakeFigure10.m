function dfdMakeFigure10()
%% Function to reproduce Figure 10 (Spatialmap) for pre vs post denoising broadband responses
%
% dfdMakeFigure10()
%
% Eline Kupers, Helena X. Wang, Kaoru Amano, Kendrick N. Kay, David J.
% Heeger, Jonathan Winawer. (YEAR) Broadband spectral responses in visual
% cortex revealed by a new MEG denoising algorithm.
% (JOURNAL. VOLUME. ISSUE. DOI.).
%
% This figure will show an interpolated spatial map of the SNR values in
% each channel for the broadband signals for each subject before and after
% using the denoising algorithm.
%
% This function assumes that data is downloaded with the DFDdownloaddata
% function.

%% Choices to make:
whichSubjects        = 1:8;
dataDir              = fullfile(dfdRootPath, 'analysis', 'data');   % Where to save data?
figureDir            = fullfile(dfdRootPath, 'analysis', 'figures');% Where to save images?
saveFigures          = true;   % Save figures in the figure folder?


%% Plot spatial map figures: right minus left stimulation for broadband component after denoising.
figure('position',[1,600,2000,300], 'Name', 'Figure 10', 'NumberTitle', 'off');

for k = 1:length(whichSubjects)
    
    whichSubject = whichSubjects(k);
    fprintf(' Load subject %d \n', whichSubject);
    data = prepareData(dataDir,whichSubject,10);    
    
    
    % Before
    subplot(2,8,k)
    
    ab_snr1 = getsignalnoise(data.results.origmodel(1),  [1 0 0], 'SNR',data.badChannels);
    ab_snr1 = to157chan(ab_snr1,~data.badChannels,'nans');
    [~,ch] = megPlotMap(ab_snr1,[-8,8],gcf,'bipolar',sprintf('S %d', k));
    set(ch,'YTick',[-8,-4,0,4,8]);
    makeprettyaxes(ch,9,9);
    
    % After
    subplot(2,8,length(whichSubjects)+k)
    ab_snr2 = getsignalnoise(data.results.finalmodel(1),  [1 0 0], 'SNR',data.badChannels);
    ab_snr2 = to157chan(ab_snr2,~data.badChannels,'nans');
    [~,ch] = megPlotMap(ab_snr2,[-8,8],gcf,'bipolar',sprintf('S %d', k));
    set(ch,'YTick',[-8,-4,0,4,8]);
    makeprettyaxes(ch,9,9);
end

if saveFigures
    hgexport(gcf,fullfile(figureDir,'figure10ab_BeforeAfter'));
%     figurewrite(fullfile(figureDir,'figure10ab_BeforeAfter'),[],0,'.',1);
end



