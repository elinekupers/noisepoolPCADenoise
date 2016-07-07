%% Calculate amplitude difference


% 1. Load results
% 2. Get top 10
% 3. Calculate signal and noise beta values
dataDir = fullfile(dfdRootPath, 'analysis','data');
allBBResults = NaN(8,2,1);
allSLResults = NaN(8,2,1);

% Define frequencies for broadband
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


for whichSubject = 1:8;
    
    fprintf('\nLoading data from subject %d of 8', whichSubject);
    [data,design] = prepareData(dataDir,whichSubject,4);
    
    % Get top ten broadband channels per subject
    load(sprintf(fullfile(dataDir, 's%02d_denoisedData_rm1epoch_bb.mat'),whichSubject));    
    topChanBB = getTop10(results);
    
    % Get top ten stimulus channels per subject
    load(sprintf(fullfile(dataDir, 's%02d_denoisedData_rm1epoch_sl.mat'),whichSubject));    
    topChanSL = getTop10(results);

    % Define conditions: Full, right, left, off
    condEpochs1 = {design{1}(:,1)==1, design{1}(:,2)==1, design{1}(:,3)==1, all(design{1}==0,2)};
    condEpochs2 = {design{2}(:,1)==1, design{2}(:,2)==1, design{2}(:,3)==1, all(design{2}==0,2)};
    
    condEpochs = {condEpochs1 condEpochs2};
    
    for dd = 1:2 % for either pre and post denoising
        
        % get power for full and blank epochs, at the specified frequencies
        this_data_full  = data{dd}(:,:,condEpochs{dd}{1});
        this_data_blank = data{dd}(:,:,condEpochs{dd}{4});
        
        bbFull  = getbroadband(this_data_full,keep_frequencies,fs);
        bbBlank = getbroadband(this_data_blank,keep_frequencies,fs);
        
        slFull  = getstimlocked(this_data_full,13).^2;  % square to get units of power
        slBlank = getstimlocked(this_data_blank,13).^2; % square to get units of power
        
        % Get mean across epochs for the top 10 channels
        %   Question 1: should we take exp(log(mean)) rather than just mean?
        %   Question 2: should we use 10 channels or fewer?
        meanFullbb  = mean(bbFull(:,topChanBB));
        meanBlankbb = mean(bbBlank(:,topChanBB));
        meanFullsl  = mean(slFull(:,topChanSL));
        meanBlanksl = mean(slBlank(:,topChanSL));
        
        
        allBBResults(whichSubject,dd,:) = mean((meanFullbb - meanBlankbb)./meanBlankbb);
        
        allSLResults(whichSubject,dd,:) = mean((meanFullsl - meanBlanksl)./meanBlanksl);
        
    end
    
end

%%
fprintf('\n');
fprintf('Mean change in broadband power, undenoised data: %4.1f%% ±  %4.1f\n', ...
    100*mean(allBBResults(:,1)), 100*std(allBBResults(:,1))/sqrt(8));

fprintf('Mean change in broadband power, denoised data: %4.1f%% ±  %4.1f\n', ...
    100*mean(allBBResults(:,2)), 100*std(allBBResults(:,2))/sqrt(8));

fprintf('Mean change in stimulus locked amplitude, undenoised data: %4.1f%% ±  %4.1f\n', ...
    100*mean(allSLResults(:,1)), 100*std(allSLResults(:,1))/sqrt(8));

% fprintf('Mean change in stimulus locked amplitude, denoised data: %4.1f%% ±  %4.1f\n', ...
%     100*mean(allSLResults(:,2)), 100*std(allSLResults(:,2))/sqrt(8));

% 
% % ECoG-BB: 2.9 -fold, 290% 0.46 log10 units
% %  MEG-BB: 0.09-fold,   9% 
% 
% 
% % ECoG BB: [1 4]
% 
% 
% % Summarize difference between 2 numbers by fold, percent difference, log10
% % units
% fold     = @(x) diff(x)./x(1,:);
% pct      = @(x) diff(x)./x(1,:)*100;
% logdelta = @(x) log10(x(2,:)./x(1,:));
% 
% betas = [1 1 1 1; 1.09 3.9 5.1 22.5];
% 
% disp(fold(betas))
% disp(pct(betas))
% disp(logdelta(betas))


