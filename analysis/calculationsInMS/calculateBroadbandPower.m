%% calculateBroadbandandSLResponses

% This is a script to calculate the relative broadband power and stimulus
% locked amplitude relative to blank periods in individual subjects, and
% the mean across subjects.

% General flow of this script:
% 0. Predefine which frequencies to use and arrays to store data
% 1. Load results
% 2. Get top 5 channels per subject
% 3. Calculate signal and noise beta values

% AUTHORS. YEAR. TITLE. JOURNAL.

%% Step 0. Predefine which frequencies to use and matrices to store data

numsubjects = 8;

topn = 5; % how many channels for each subject?

% Define data directory and predefine arrays
dataDir = fullfile(dfdRootPath, 'analysis','data');
allBBResults = NaN(numsubjects,2,1);
allSLResults = NaN(numsubjects,2,1);

% Define frequencies to compute the broadband power
fs           = 1000;         % Sample rate
f            = 0:150;        % Limit frequencies to [0 150] Hz
slFreq       = 12;           % Stimulus-locked frequency
tol          = 1.5;          % Exclude frequencies within +/- tol of sl_freq
slDrop       = f(mod(f, slFreq) <= tol | mod(f, slFreq) > slFreq - tol);

% Exclude all frequencies below 60 Hz when computing broadband power
lfDrop       = f(f<60);

% Define the frequenies and indices into the frequencies used to compute
% broadband power
[~, abIndex] = setdiff(f, [slDrop lfDrop]);

% Create function handles for the frequencies that we use
keepFrequencies    = @(x) x(abIndex);
bbFrequencies      = f(abIndex);

% Note the geometric mean of the broadband frequencies. The broadband
% metric is equivalent to fitting a line in log power/log frequency space,
% and evaluating this line at the geometric mean of the broadband
% frequencies.
bb_eval = geomean(f(abIndex));

for whichSubject = 1:8
    
    %% Step 1. Load results
    fprintf('\nLoading data from subject %d of 8', whichSubject);
    [data,design] = prepareData(dataDir,whichSubject,4); % Figure 4 contains the same data, therefore we prepare data from this figure
    
    %% Step 2a: Get top five broadband channels per subject
    load(sprintf(fullfile(dataDir, 's%02d_denoisedData_bb.mat'),whichSubject));    
    topChanBB = getTop10(results, [], topn);
    
    % Step 2b: Get top five stimulus channels per subject
    load(sprintf(fullfile(dataDir, 's%02d_denoisedData_sl.mat'),whichSubject));    
    topChanSL = getTop10(results, [], topn);

    % Define conditions: Full, right, left, off
    condEpochs1 = {design{1}(:,1)==1, design{1}(:,2)==1, design{1}(:,3)==1, all(design{1}==0,2)};
    condEpochs2 = {design{2}(:,1)==1, design{2}(:,2)==1, design{2}(:,3)==1, all(design{2}==0,2)};
    
    condEpochs = {condEpochs1 condEpochs2};
    
    for dd = 1:2 % for either 1:pre and 2:post denoising
        

        %% Step 3: Compute log power for full (1) and blank (4) epochs, at the specified frequencies
        this_data_full  = data{dd}(:,:,condEpochs{dd}{1});
        this_data_blank = data{dd}(:,:,condEpochs{dd}{4});
        
        bbFull  = log10(getbroadband(this_data_full,keepFrequencies,fs)); % Broadband data is already in units of power
        bbBlank = log10(getbroadband(this_data_blank,keepFrequencies,fs)); % Broadband data is already in units of power
        
        slFull  = log10(getstimlocked(this_data_full,13).^2);  % Square to get units of power
        slBlank = log10(getstimlocked(this_data_blank,13).^2); % Square to get units of power
        
        % Mean across epochs for each of the top 5 channels
        meanFullbb  = mean(bbFull(:,topChanBB));
        meanBlankbb = mean(bbBlank(:,topChanBB));
        meanFullsl  = mean(slFull(:,topChanSL));
        meanBlanksl = mean(slBlank(:,topChanSL));
       
        % Full field minus blank, averaged across top 5 channels
        allBBResults(whichSubject,dd,:) = mean(meanFullbb - meanBlankbb);
        
        % Add individual subjects stimulus locked results to all results
        allSLResults(whichSubject,dd,:) = mean(meanFullsl - meanBlanksl);
        
    end
    
end

%%
fprintf('\n');

meg.bb_pre.mn = 100*(mean(10.^allBBResults(:,1))-1);
meg.bb_pre.se = 100*std(10.^allBBResults(:,1))/sqrt(numsubjects);

meg.bb_post.mn = 100*(mean(10.^allBBResults(:,2))-1);
meg.bb_post.se = 100*std(10.^allBBResults(:,2))/sqrt(numsubjects);

meg.sl.mn = 100*(mean(10.^allSLResults(:,1))-1);
meg.sl.se = 100*std(10.^allSLResults(:,1))/sqrt(numsubjects);

fprintf('Mean change in broadband power, undenoised data: %4.1f%% ± %3.1f%%\n', ...
    meg.bb_pre.mn, meg.bb_pre.se);


fprintf('Mean change in broadband power, denoised data: %4.1f%% ± %3.1f%%\n', ...
    meg.bb_post.mn, meg.bb_post.se);

fprintf('Mean change in stimulus locked amplitude, undenoised data: %4.1f%% ±  %3.1f%%\n', ...
   meg.sl.mn, meg.sl.se);

%% Do same calculations for ECoG
ecog = load('ECoG_OnOffData');

%% broadband

ecog.bbFull  = log10(ecog.data.mnOn(:,1));
ecog.bbBlank = log10(ecog.data.mnOff(:,1));
ecog.allBBResults = ecog.bbFull - ecog.bbBlank;

ecog.slFull  = log10(ecog.data.mnOn(:,2));
ecog.slBlank = log10(ecog.data.mnOff(:,2));
ecog.allSLResults = ecog.slFull - ecog.slBlank;

numecog = length(ecog.bbFull);


ecog.bb.mn = 100*(mean(10.^ecog.allBBResults)-1);
ecog.bb.se = 100*std(10.^ecog.allBBResults)/sqrt(numecog);

ecog.sl.mn = 100*(mean(10.^ecog.allSLResults)-1);
ecog.sl.se = 100*std(10.^ecog.allSLResults)/sqrt(numecog);



fprintf('Mean change in broadband power, ECoG data: %4.1f%% ± %3.1f%%\n', ...
    ecog.bb.mn, ecog.bb.se);

fprintf('Mean change in stimulus locked amplitude, ECoG data: %4.1f%% ± %3.1f%%\n', ...
    ecog.sl.mn, ecog.sl.se);

%% Ratios

 
MEG_sl_to_bb  = meg.sl.mn / meg.bb_post.mn
ECoG_sl_to_bb = ecog.sl.mn / ecog.bb.mn

SL_ECoG_to_MEG = ecog.sl.mn / meg.sl.mn
BB_ECoG_to_MEG =  ecog.bb.mn / meg.bb_post.mn


