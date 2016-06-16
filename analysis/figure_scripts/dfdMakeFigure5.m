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
whichSubject    = 1;        % Subject 1 is the example subject.
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
figure('position',[1,600,1400,800]); set(gcf, 'Name', 'Figure 5, Example subject');
condNames = {'Stim Full','Stim Left','Stim Right'};
for icond = 1:3
    % get stimulus-locked snr
    sl_snr1 = getsignalnoise(sl.results.origmodel(1),icond, 'SNR',sl.badChannels);
    % get broadband snr for before and after denoising
    ab_snr1 = getsignalnoise(bb.results.origmodel(1),  icond, 'SNR',bb.badChannels);
    
    % get difference SNR
    sl_snr1_L = getsignalnoise(sl.results.origmodel(1), 2, 'SNR');
    sl_snr1_R = getsignalnoise(sl.results.origmodel(1), 3, 'SNR');

    ab_snr1_L = getsignalnoise(bb.results.origmodel(1),  2, 'SNR');
    ab_snr1_R = getsignalnoise(bb.results.origmodel(1),  3, 'SNR');


    if whichSubject < 9; % NeuroMag360 data is already converted when combining the channels
        % convert back into 157-channel space  
        sl_snr1 = to157chan(sl_snr1,~sl.badChannels,'nans');
        ab_snr1 = to157chan(ab_snr1,~bb.badChannels,'nans');
        
        sl_snr1a_LmnR = to157chan(sl_snr1_L,~bb.badChannels,'nans') - to157chan(sl_snr1_R,~bb.badChannels,'nans');
        ab_snr1a_LmnR = to157chan(ab_snr1_L,~bb.badChannels,'nans') - to157chan(ab_snr1_R,~bb.badChannels,'nans');
    end
    
    % Threshold
    ab_snr1(abs(ab_snr1) < threshold) = 0;
    sl_snr1(abs(sl_snr1) < threshold) = 0;
    
    ab_snr1a_LmnR(abs(ab_snr1a_LmnR) < threshold) = 0;
    sl_snr1a_LmnR(abs(sl_snr1a_LmnR) < threshold) = 0;
    
   %CHECK IF THIS IS THE SAME 
%     ab_beta2 = results.whichmodel.beta(1,:,:); %right
%     ab_diff = ab_beta2 - ab_beta1;
%     
%     diff_med = nanmedian(squeeze(ab_diff),2);
%     diff_se = nanstd(squeeze(ab_diff),[],2);
%         
%     ab_snr_diff = to157chan((diff_med./diff_se)',~data{1}.badChannels,'nans');
%     
%     subplot(2,4,k);  

    % Set colormap limits
%     max_val = max(abs([ab_snr1a_LmnR, ab_snr2a_LmnR]));
%     clims_ab = [-1,1].*[max_val,max_val];
    clims_sl = [-25.6723,25.6723];
    clims_ab = [-8.4445, 8.4445];
    clims_ab_diff = [-5.5363, 5.5363];

    % plot spatial maps
    subplot(4,2,(icond-1)*2+1)
    [~,ch] = megPlotMap(sl_snr1,clims_sl,gcf,'bipolar',sprintf('%s : Stimulus Locked Original', condNames{icond}),data_hdr,cfg);
    makeprettyaxes(ch,9,9);
    set(ch,'YTick',[-20,-10,0,10,20]);
    title(sprintf('SL no DN %s', condNames{icond}))
    
    subplot(4,2,(icond-1)*2+2)
    [~,ch] = megPlotMap(ab_snr1,clims_ab,gcf,'bipolar',sprintf('%s Original', condNames{icond}),data_hdr,cfg);
    makeprettyaxes(ch,9,9);
    set(ch,'YTick',[-8,-4,0,4,8]);
    title(sprintf('Broadband Pre %s', condNames{icond}))
    
    % plot difference spatial maps
    subplot(4,2,7)
    [~,ch] = megPlotMap(sl_snr1a_LmnR,clims_sl,gcf,'bipolar',sprintf('%s : Stimulus Locked Original', 'Left minus Right'));
    makeprettyaxes(ch,9,9);
    set(ch,'YTick',[-20,-10,0,10,20]);
    title(sprintf('SL no DN %s', 'Left minus Right'))
    
    subplot(4,2,8)
    [~,ch] = megPlotMap(ab_snr1a_LmnR,clims_ab_diff,gcf,'bipolar',sprintf('%s Original', 'Left minus Right'));
    makeprettyaxes(ch,9,9);
    set(ch,'YTick',[-5,-2.5,0,2.5,5]);
    title(sprintf('Broadband Pre %s', 'Left minus Right'))
    
end

if saveFigures
    printnice(gcf, 0, figureDir, sprintf('figure5_examplesubject%d_bipolar_thresh%d_interpolated',whichSubject, threshold));
end

%% Now call dfdMakeFigure5AcrossSubjects
dfdMakeFigure5AcrossSubjects();


