%% Calculate amplitude difference


% 1. Load results
% 2. Get top 10
% 3. Calculate signal and noise beta values
dataDir = fullfile(dfdRootPath, 'analysis','data');
allBBResults = NaN(8,2,1);
allSLResults = NaN(8,2,1);


for whichSubject = 1:8;
    
    [data,design,exampleIndex,exampleChannel] = prepareData(dataDir,whichSubject,4);
    
    load(sprintf(fullfile(dataDir, 's%02d_denoisedData_rm1epoch_bb.mat'),whichSubject));
    
    % Get top ten channels per subject
    topChan = getTop10(results);
    
    % Define conditions: Full, right, left, off
    condEpochs1 = {design{1}(:,1)==1, design{1}(:,2)==1, design{1}(:,3)==1, all(design{1}==0,2)};
    condEpochs2 = {design{2}(:,1)==1, design{2}(:,2)==1, design{2}(:,3)==1, all(design{2}==0,2)};

    condEpochs = {condEpochs1 condEpochs2};

    
    %% Define frequencies for broadband
    fs          = 1000;
    f           = 0:150;   % limit frequencies to [0 150] Hz
    sl_freq     = 12;      % Stimulus-locked frequency
    sl_freq_i   = sl_freq + 1;
    tol         = 1.5;     % exclude frequencies within +/- tol of sl_freq
    sl_drop     = f(mod(f, sl_freq) <= tol | mod(f, sl_freq) > sl_freq - tol);
    
    % Exclude all frequencies that are close to a multiple of the
    % line noise frequency (60 Hz)
    ln_drop     = f(mod(f, 60) <= tol | mod(f, 60) > 60 - tol);
    
    % Exclude all frequencies below 60 Hz when computing broadband power
    lf_drop     = f(f<60);
    
    % Define the frequenies and indices into the frequencies used to compute
    % broadband power
    [~, ab_i]   = setdiff(f, [sl_drop ln_drop lf_drop]);
    
    keep_frequencies    = @(x) x(ab_i);
    bb_frequencies      = f(ab_i);
    
    for dd = 1:2 % for either pre and post denoising
        
        % get power for full and blank epochs, at the specified frequencies
    	this_data_full = data{dd}(:,:,condEpochs{dd}{1});
        this_data_blank = data{dd}(:,:,condEpochs{dd}{4});
        
        bbFull = getbroadband(this_data_full,keep_frequencies,fs);
        bbBlank = getbroadband(this_data_blank,keep_frequencies,fs);
        
        slFull = getstimlocked(this_data_full,13);
        slBlank = getstimlocked(this_data_blank,13);
        
           % Get mean across top 10 Across epochs
         allBBResults(whichSubject,dd,:) = mean(mean(bbFull(:,topChan)) -mean(bbBlank(:,topChan)));

         allSLResults(whichSubject,dd,:) = mean(mean(slFull(:,topChan)) -mean(slBlank(:,topChan)));

    end
    
end

disp(mean(allBBResults(:,1)));
disp(mean(allBBResults(:,2)));

disp(std(allBBResults(:,1))/sqrt(8));
disp(std(allBBResults(:,2))/sqrt(8));

disp(mean(allSLResults(:,1)));
disp(mean(allSLResults(:,2)));

disp(std(allSLResults(:,1))/sqrt(8));
disp(std(allSLResults(:,2))/sqrt(8));


