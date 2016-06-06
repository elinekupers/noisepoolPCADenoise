function dfdMakeFigure10()
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
figure('position',[1,600,1400,800]);
condNames = {'Stim Full','Stim Left','Stim Right'};
for icond = 1:3
    % get broadband snr for before and after denoising
    ab_snr1 = getsignalnoise(bb.results.origmodel(1),  icond, 'SNR',bb.badChannels);
    ab_snr2 = getsignalnoise(bb.results.finalmodel(1), icond, 'SNR',bb.badChannels);
    
    ab_snr1_L = getsignalnoise(bb.results.origmodel(1),  2, 'SNR');
    ab_snr1_R = getsignalnoise(bb.results.origmodel(1),  3, 'SNR');

    ab_snr2_L = getsignalnoise(bb.results.finalmodel(1), 2, 'SNR');
    ab_snr2_R = getsignalnoise(bb.results.finalmodel(1), 3, 'SNR');

    if whichSubject < 9; % NeuroMag360 data is already converted when combining the channels
        % convert back into 157-channel space  
        ab_snr1 = to157chan(ab_snr1,~bb.badChannels,'nans');
        ab_snr2 = to157chan(ab_snr2,~bb.badChannels,'nans');

        ab_snr1a_LmnR = to157chan(ab_snr1_L,~bb.badChannels,'nans') - to157chan(ab_snr1_R,~bb.badChannels,'nans');
        ab_snr2a_LmnR = to157chan(ab_snr2_L,~bb.badChannels,'nans') - to157chan(ab_snr2_R,~bb.badChannels,'nans');
    end
    
    % Threshold
    ab_snr1(abs(ab_snr1) < threshold) = 0;
    ab_snr2(abs(ab_snr2) < threshold) = 0;
    
    ab_snr1a_LmnR(abs(ab_snr1a_LmnR) < threshold) = 0;
    ab_snr2a_LmnR(abs(ab_snr2a_LmnR) < threshold) = 0;
    
    % Set colormap limits
%     max_val = max(abs([ab_snr1a_LmnR, ab_snr2a_LmnR]));
%     clims_ab = [-1,1].*[max_val,max_val];
    clims_sl = [-25.6723,25.6723];
    clims_ab = [-8.4445, 8.4445];     

    subplot(4,2,(icond-1)*2+1)
    [~,ch] = megPlotMap(ab_snr1,clims_ab,gcf,'bipolar',sprintf('%s Original', condNames{icond}),data_hdr,cfg);
    makeprettyaxes(ch,9,9);
    title(sprintf('Broadband Pre %s', condNames{icond}))
    
    subplot(4,2,(icond-1)*2+2)
    [~,ch] = megPlotMap(ab_snr2,clims_ab,gcf,'bipolar',sprintf('%s : Denoised PC %d',condNames{icond}, bb.results.pcnum(1)),data_hdr,cfg);
    makeprettyaxes(ch,9,9);
    title(sprintf('Broadband Post %s', condNames{icond}))
    
    % plot difference spatial maps
    subplot(4,2,(icond-1)*2+3)
    [~,ch] = megPlotMap(ab_snr1a_LmnR,clims_ab,gcf,'bipolar',sprintf('%s Original', 'Left minus Right'));
    makeprettyaxes(ch,9,9);
    set(ch,'YTick',[-5,-2.5,0,2.5,5]);
    title(sprintf('Broadband Pre %s', 'Left minus Right'))
    
    subplot(4,2,(icond-1)*2+4)
    [~,ch] = megPlotMap(ab_snr2a_LmnR,clims_ab,gcf,'bipolar',sprintf('%s : Denoised PC %d', 'Left minus Right', bb.results.pcnum(1)));
    makeprettyaxes(ch,9,9);
    set(ch,'YTick',[-5,-2.5,0,2.5,5]);
    title(sprintf('Broadband Post %s', 'Left minus Right'))
    
end

if saveFigures
    printnice(gcf, 0, figureDir, sprintf('figure5_examplesubject%d_bipolar_thresh%d_interpolated',whichSubject, threshold));
end

%% Now do the same but then across subjects
dfdMakeFigure10AcrossSubjects()

