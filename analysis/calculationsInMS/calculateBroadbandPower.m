%% calculateBroadbandPower

% This is a script to calculate the relative broadband power compared to blank
% periods in individual subjects spectra and the mean across subjects. 

% General flow of this script:
% 0. Predefine which frequencies to use and arrays to store data
% 1. Load results
% 2. Get top 5 channels per subject
% 3. Calculate signal and noise beta values

% AUTHORS. YEAR. TITLE. JOURNAL.

%% Step 0. Predefine which frequencies to use and matrices to store data

% Define data directory and predefine arrays
dataDir = fullfile(dfdRootPath, 'analysis','data');
allBBResults = NaN(8,2,1);
allSLResults = NaN(8,2,1);

% Define frequencies to compute the broadband power
fs           = 1000;         % Sample rate
f            = 0:150;        % Limit frequencies to [0 150] Hz
slFreq      = 12;           % Stimulus-locked frequency
% slFreqIndex   = slFreq + 1;  % Stimulus-locked index
tol          = 1.5;          % Exclude frequencies within +/- tol of sl_freq
slDrop       = f(mod(f, slFreq) <= tol | mod(f, slFreq) > slFreq - tol);

% Exclude all frequencies that are close to a multiple of the
% line noise frequency (60 Hz)
lnDrop       = f(mod(f, 60) <= tol | mod(f, 60) > 60 - tol);

% Exclude all frequencies below 60 Hz when computing broadband power
lfDrop       = f(f<60);

% Define the frequenies and indices into the frequencies used to compute
% broadband power
[~, abIndex] = setdiff(f, [slDrop lnDrop lfDrop]);

% Create function handles for the frequencies that we use
keepFrequencies    = @(x) x(abIndex);
bbFrequencies      = f(abIndex);


for whichSubject = 1:8;
    
    %% Step 1. Load results
    fprintf('\nLoading data from subject %d of 8', whichSubject);
    [data,design] = prepareData(dataDir,whichSubject,4); % Figure 4 contains the same data, therefore we prepare data from this figure
    
    %% Step 2a: Get top five broadband channels per subject
    load(sprintf(fullfile(dataDir, 's%02d_denoisedData_rm1epoch_bb.mat'),whichSubject));    
    topChanBB = getTop10(results, [], 5);
    
    % Step 2b: Get top five stimulus channels per subject
    load(sprintf(fullfile(dataDir, 's%02d_denoisedData_rm1epoch_sl.mat'),whichSubject));    
    topChanSL = getTop10(results, [], 5);

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
        
        % Get mean across epochs for the top 10 channels
        meanFullbb  = mean(bbFull(:,topChanBB));
        meanBlankbb = mean(bbBlank(:,topChanBB));
        meanFullsl  = mean(slFull(:,topChanSL));
        meanBlanksl = mean(slBlank(:,topChanSL));
       
        % Add individual subjects broadband results to all results
        allBBResults(whichSubject,dd,:) = mean((meanFullbb - meanBlankbb));
        
        % Add individual subjects stimulus locked results to all results
        allSLResults(whichSubject,dd,:) = mean((meanFullsl - meanBlanksl));
        
    end
    
end

%%
fprintf('\n');

fprintf('Mean change in broadband power, undenoised data: %4.4f ±  %4.4f\n', ...
    mean(allBBResults(:,1)), std(allBBResults(:,1))/sqrt(8));

fprintf('Mean change in broadband power, denoised data: %4.4f ±  %4.4f\n', ...
    mean(allBBResults(:,2)), std(allBBResults(:,2))/sqrt(8));

fprintf('Mean change in stimulus locked amplitude, undenoised data: %4.4f ±  %4.4f\n', ...
    mean(allSLResults(:,1)), std(allSLResults(:,1))/sqrt(8));

