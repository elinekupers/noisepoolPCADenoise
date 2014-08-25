%% Define paths and data sets 
inputDataDir = '/Volumes/HelenaBackup/denoisesuite/tmpmeg/';
fitDataStr = 'b2fr_hpf2_fit10p1k';
%fitDataStr = 'b2frSL_fit10p1k';
whichfun   = 1;
figuredir = 'manuscript_figs/figure_spatialmap';

%% example subject spatial map, SL: L,R, Full, Broadband: L,R, Full - Fig.8
% Define example session 
sessionNum = 3;
sessionDir = megGetDataPaths(sessionNum);
% stimulus locked glm results 
slresults = load(fullfile(inputDataDir,sprintf('%sb2fSL_fitfull75',sessionDir)));
slresults = slresults.results;
% broadband glm results 
thisfile = fullfile(inputDataDir,sprintf('%s%s',sessionDir,fitDataStr));
disp(thisfile); load(thisfile,'results');
% load 'badChannels' matrix
datafile = fullfile(inputDataDir,'inputdata',sprintf('%sb2',sessionDir));
disp(datafile); load(datafile,'badChannels');

%% Plot the subject 
figure('position',[1,600,1400,800]);
for icond = 1:3
    % get stimulus-locked snr 
    sl_snr1 = getsignalnoise(slresults.origmodel(1),icond, 'SNR');
    %clims_sl = [0, max(sl_snr1)];
    clims_sl = [0,25.6723];
    % get broadband snr for before and after denoising 
    ab_snr1 = getsignalnoise(results.origmodel(whichfun),  icond, 'SNR');
    ab_snr2 = getsignalnoise(results.finalmodel(whichfun), icond, 'SNR');
    clims_ab = [0, max([ab_snr1, 12.4445])];
    %clims_ab = [0, max([ab_snr1, ab_snr2])];
    
    % convert back into 157-channel space 
    ab_snr1a = to157chan(ab_snr1,~badChannels,'nans');
    ab_snr2a = to157chan(ab_snr2,~badChannels,'nans');
    sl_snr1a = to157chan(sl_snr1,~badChannels,'nans');
    
    % plot spatial maps 
    subplot(3,3,(icond-1)*3+1)
    [~,ch] = megPlotMap(sl_snr1a,clims_sl,gcf,'jet','Stimulus Locked Original');
    makeprettyaxes(gca,9,9);
    makeprettyaxes(ch,9,9);
    
    subplot(3,3,(icond-1)*3+2)
    [~,ch] = megPlotMap(ab_snr1a,clims_ab,gcf,'jet','Original');
    makeprettyaxes(gca,9,9);
    makeprettyaxes(ch,9,9);
    
    subplot(3,3,(icond-1)*3+3)
    [~,ch] = megPlotMap(ab_snr2a,clims_ab,gcf,'jet',sprintf('Denoised PC %d',results.pcnum(whichfun)));
    makeprettyaxes(gca,9,9);
    makeprettyaxes(ch,9,9);
end

%figurewrite(fullfile(figuredir,'figure_s3'),[],0,'.',1);

%% Plot all subjects - full condition, post minus pre - Fig. 9

figure('position',[1,600,1400,800]);
sessionNums = [11,12, 3:6, 9:10];%[1:6,9,10];
icond = 1;

for k = 1:length(sessionNums)
    sessionDir = megGetDataPaths(sessionNums(k));
    thisfile = fullfile(inputDataDir,sprintf('%s%s',sessionDir,fitDataStr));
    disp(thisfile); load(thisfile,'results');
    datafile = fullfile(inputDataDir,'inputdata',sprintf('%sb2',sessionDir));
    disp(datafile); load(datafile,'badChannels');
    
    subplot(2,4,k);
    ab_snr1 = getsignalnoise(results.origmodel(whichfun),  icond, 'SNR');
    ab_snr2 = getsignalnoise(results.finalmodel(whichfun), icond, 'SNR');
    ab_snr_diff = to157chan(ab_snr2-ab_snr1,~badChannels,'nans');
    
    [~,ch] = megPlotMap(ab_snr_diff,[-10,10],gcf,jmaColors('coolhotcortex'));
    makeprettyaxes(gca,9,9);
    set(ch,'ytick',-10:5:10); 
    makeprettyaxes(ch,9,9); 
    %keyboard
end
%figurewrite(fullfile(figuredir,'figure_bbdiffall'),[],0,'.',1);

%% All subjects - right minus left - Fig. 9b

figure('position',[1,600,1400,800]);
%whichmodel = 'origmodel';
whichmodel = 'finalmodel';
for k = 1:length(sessionNums)
    sessionDir = megGetDataPaths(sessionNums(k));
    thisfile = fullfile(inputDataDir,sprintf('%s%s',sessionDir,fitDataStr));
    disp(thisfile); load(thisfile,'results');
    datafile = fullfile(inputDataDir,'inputdata',sprintf('%sb2',sessionDir));
    disp(datafile); load(datafile,'badChannels');
    
    subplot(2,4,k);
%     ab_snr1 = getsignalnoise(results.origmodel(whichfun), 2, 'SNR');
%     ab_snr2 = getsignalnoise(results.finalmodel(whichfun),2, 'SNR');
%     ab_snr_diff1 = to157chan(ab_snr2-ab_snr1,~badChannels,'nans');
%     ab_snr1 = getsignalnoise(results.origmodel(whichfun), 3, 'SNR');
%     ab_snr2 = getsignalnoise(results.finalmodel(whichfun),3, 'SNR');
%     ab_snr_diff2 = to157chan(ab_snr2-ab_snr1,~badChannels,'nans');
%     ab_snr_diff = ab_snr_diff1-ab_snr_diff2;

    ab_snr1 = getsignalnoise(results.(whichmodel)(whichfun), 2, 'SNR');
    ab_snr2 = getsignalnoise(results.(whichmodel)(whichfun), 3, 'SNR');
    ab_snr_diff = to157chan(ab_snr2-ab_snr1,~badChannels,'nans');
    
    [~,ch] = megPlotMap(ab_snr_diff,[-5,5],gcf,jmaColors('coolhotcortex'));
    makeprettyaxes(gca,9,9);
    %set(ch,'ytick',-20:10:20);
    set(ch,'ytick',-5:5:5);
    makeprettyaxes(ch,9,9); 
    %keyboard
end
%figurewrite(fullfile(figuredir,'figure_bbRightMLeft_after'),[],0,'.',1);