%% Plot multiple denoised sessions of one subject
%
% We do this in order to check how bootstrapping influences our final
% results. We will do this for session with name 8, so session number 5.

%% Choices to make:

nr_of_repeats      = 5;

sessionNums        = 5;   % You need all subjects for this figure 
                            % Choose a particular number if you would like
                            % a specific subject
plotBb             = true;  % True = Broadband, false = Stim locked                            

inputDataDir       = fullfile(DFDrootpath, 'data'); % Note: New data matrices will 
                                                    % also get stored in the same folder.
whichFun           = 1;     % ToDo: figure out what this means..

conditionNumbers   = 1:6;   % Choose 1:6 to get all conditions: Full, 
                            % left, right (for 'on' 1,2,3 and 'off' 4,5,6)
                            % conditions. 
                            % If you like only a certain conditions e.g.
                            % full on and off define this variable as [1,4]

sensorDataStr      = 'b2';  % Use second version of defining data:
                            % data created by allowing fewer 0 epochs, 
                            % using the new data format (such that design 
                            % matrix and condition order are as they occurred 
                            % in the experiment)
                            
freqDefinition     = 'f';   % f  - uses new definition of ab_i for freq 
                            % (includes more points; excludes 1 pt either side rather than 3)

noisePoolSelection = 'r';   % r - noise pool selection (and pc cutoff, if 
                            % selected by algorithm) by SNR rather than by R2

doHPF              = 'hpf2';% hpf2 - high pass filtered with a sharp cutoff
                            % at 62 Hz and without stimulus harmonics                        
                            
nrPCs              = 'fit10'; % fit10 - jump to 10 PCs as cutoff,
                            % In Process: (fitfull10 - don't jump to 10, but denoise all
                            %             channels in between 0 and 10
                            %             PCs.)
                            % fitfull75 - use all channels in noisepool   

xBoots             = 'p1k'; % p1k - bootstrapped 1000x rather than 100x (p100)

fitDataStrBB       = [sensorDataStr freqDefinition noisePoolSelection '_' ...
                        doHPF '_' nrPCs xBoots];
                    
fitDataStrSL       = [sensorDataStr freqDefinition noisePoolSelection 'SL_' ...
                        nrPCs xBoots]; 

saveData         = true;    % Separate matfiles are saved, in order to 
                            % speed up the script if you only want to plot.
saveEpochGroup   = false;   % Epochs can be grouped in a certain order, 
                            % you can save this if you like.
                            
figureDir        = fullfile(DFDrootpath, 'figures');

saveFigures      = true;   % Save figures in the figure folder?



%% Check whether we got our preprocessed data matrices

% % Get session name and top directory
% % This can be modified so that top directory points somewhere else
% dataset = DFDgetdatapaths(sessionNums,conditionNumbers,inputDataDir);
% if ~exist(fullfile(inputDataDir, 'savedProcData', [dataset fitDataStrBB '.mat']),'file'); end
% % if ~exist(fullfile(inputDataDir, 'savedProcData', [dataset fitDataStrSL '.mat']),'file'); end
% 
% [sensorData, design, badChannels, conditionNames, okEpochs] = ...
%     DFDpreload(sessionNums, sensorDataStr, saveData, saveEpochGroup, inputDataDir, conditionNumbers);
% 
% %% Denoise BB
% 
% doHpc            = true; % High pass filter data 
% evalfunToCompute = {'bb'}; % Broadband
% saveDenoiseTs    = true; % You need denoised ts to make the spectrum figure.
% pcstop10         = true; % To use all PC's in the noise pool (not just number 10)
% 
% resultsBB        = DFDDenoiseWrapper(sessionNums, [], [], doHpc, [], pcstop10, ...
%                                         evalfunToCompute, [], saveDenoiseTs);

s1 = fullfile(DFDrootpath,'data','savedProcData','SSMEG_Dataset_08b2fr_hpf2_fit10p1k_1.mat');
s2 = fullfile(DFDrootpath,'data','savedProcData','SSMEG_Dataset_08b2fr_hpf2_fit10p1k_2.mat');
s3 = fullfile(DFDrootpath,'data','savedProcData','SSMEG_Dataset_08b2fr_hpf2_fit10p1k_3.mat');
s4 = fullfile(DFDrootpath,'data','savedProcData','SSMEG_Dataset_08b2fr_hpf2_fit10p1k_4.mat');

replicatedsessions = {s1; s2; s3; s4};

%% Plot it!

condColors   = [63, 121, 204; 228, 65, 69; 116,183,74]/255;

%% S, N, and SNR shown separately, before versus after denoising with 10
%% PCs. For 3 example sessions - Fig. 7A,B,C

% exampleSessions = [3,4,5]; % Helena's [5,6,9]
% for k = 1:length(exampleSessions)
%     fprintf(' session %d \n', exampleSessions(k));
%     sessionDir = DFDgetdatapaths(exampleSessions(k),conditionNumbers,inputDataDir);
    for ii = 1:length(replicatedsessions)
    % load fit file
    thisfile = replicatedsessions{ii};
    disp(thisfile); load(thisfile,'results');
    
    % top 10 channels
    pcchan = getTop10(results,whichFun);
    %pcchan = results.pcchan{whichFun};
    %pcchan = ~results.noisepool;
    %pcchanfile = fullfile(inputDataDir,sprintf('%s%s',sessionDir,'b2fr_hpf2_fitfull75'));
    %tmp = load(pcchanfile); pcchan = tmp.results.pcchan{whichFun};
    
    % signal and noise before denoising
    ab_signal1 = abs(results.origmodel(whichFun).beta_md(:,pcchan));
    ab_noise1  = results.origmodel(whichFun).beta_se(:,pcchan);
    % signal and noise after denoising 
    ab_signal2 = abs(results.finalmodel(whichFun).beta_md(:,pcchan));
    ab_noise2  = results.finalmodel(whichFun).beta_se(:,pcchan);
    % snr = signal/denoise
    ab_snr1    = ab_signal1./ab_noise1;
    ab_snr2    = ab_signal2./ab_noise2;
    
    % plot each condition as a different color 
    % signal 
    subplot(3,length(replicatedsessions),ii); cla; hold on;
    for nn = 1:3
        plot(ab_signal1(nn,:),ab_signal2(nn,:),'o','color',condColors(nn,:));
    end
    axis square;
    axismax = max([ab_signal1(:);ab_signal2(:)])*1.2;
    xlim([0,axismax]); ylim([0,axismax]); line([0,axismax],[0,axismax],'color','k');
    title(sprintf('S%d : signal', 5));
    makeprettyaxes(gca,9,9);
    
    % noise
    subplot(3,length(replicatedsessions),ii+length(replicatedsessions)); cla; hold on;
    for nn = 1:3
        plot(ab_noise1(nn,:),ab_noise2(nn,:),'o','color',condColors(nn,:));
    end
    axismax = max([ab_noise1(:); ab_noise2(:)])*1.2;
    xlim([0,axismax]); ylim([0,axismax]); line([0,axismax],[0,axismax],'color','k');
    axis square;
    title(sprintf('S%d : noise', 5));
    makeprettyaxes(gca,9,9);
    
    % snr
    subplot(3,length(replicatedsessions),ii+2*length(replicatedsessions)); cla; hold on;
    for nn = 1:3
        plot(ab_snr1(nn,:),ab_snr2(nn,:),'o','color',condColors(nn,:));
    end
    axismax = max([ab_snr1(:); ab_snr2(:)])*1.2;
    xlim([0,axismax]); ylim([0,axismax]); line([0,axismax],[0,axismax],'color','k');
    axis square;
    title(sprintf('S%d : SNR', 5));
    makeprettyaxes(gca,9,9);
    
    drawnow;
end

if saveFigures
    if plotBb
        figurewrite(fullfile(figureDir,'Figure7abc_snrexamples_Subject5compare4runs'),[],0,'.',1);
    else
        figurewrite(fullfile(figureDir,'Figure12abc_snrexamples_SL'),[],0,'.',1);
    end
end