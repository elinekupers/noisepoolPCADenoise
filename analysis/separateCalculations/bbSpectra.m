%% Calculate amplitude difference

lgmn = @(x,n) exp(mean(log(x),n));

% 1. Load results
% 2. Get top n
% 3. Calculate spectra
numSubjects  = 8;
numChannels  = 5;

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
numFrequencies      = numel(bb_frequencies);

dataDir = fullfile(dfdRootPath, 'analysis','data');
spectraFull  = cell(numSubjects,2);
spectraBlank = cell(numSubjects,2);
spectraRight = cell(numSubjects,2);
spectraLeft  = cell(numSubjects,2);

%% loop over subjects
for whichSubject = 1:numSubjects;
    
    fprintf('\nLoading data from subject %d of 8', whichSubject); drawnow()
    [data,design] = prepareData(dataDir,whichSubject,4);
    
    
    % Get top n broadband channels per subject
    load(sprintf(fullfile(dataDir, 's%02d_denoisedData_rm1epoch_bb.mat'),whichSubject));    
    topN = getTop10(results, [], numChannels);
    

    % Define conditions: Full, right, left, off
    condEpochs1 = {design{1}(:,1)==1, design{1}(:,2)==1, design{1}(:,3)==1, all(design{1}==0,2)};
    condEpochs2 = {design{2}(:,1)==1, design{2}(:,2)==1, design{2}(:,3)==1, all(design{2}==0,2)};
    
    condEpochs = {condEpochs1 condEpochs2};
    
    for dd = 1:2 % for either pre and post denoising

        pxx = (abs(fft(data{dd},[], 2))/size(data{dd},2)*2).^2;

        % get power for full and blank epochs, at the specified frequencies
        tmp = pxx(topN,ab_i,condEpochs{dd}{1});
        tmp = squeeze(lgmn(tmp,1));
        spectraFull{whichSubject, dd} = tmp;
        
        tmp = pxx(topN,ab_i,condEpochs{dd}{2});
        tmp = squeeze(lgmn(tmp,1));
        spectraRight{whichSubject, dd} = tmp;

        tmp = pxx(topN,ab_i,condEpochs{dd}{3});
        tmp = squeeze(lgmn(tmp,1));
        spectraLeft{whichSubject, dd} = tmp;
               
        tmp = pxx(topN,ab_i,condEpochs{dd}{4});
        tmp = squeeze(lgmn(tmp,1));
        spectraBlank{whichSubject, dd} = tmp;

        
    end
    
end

%% Plot 'em
 
colors = dfdGetColors(4); str = {'Before denoising' 'After denoising'};

% % Before denoising
% for dd = 1:2
%     
%     fH = figure(dd); clf, set(gcf, 'Color', 'w')
%     
%     mn = lgmn(spectraFull{dd},2);
%     se = std(spectraFull{dd},[], 2) / sqrt(size(spectraFull{dd},2));    
%     fill([bb_frequencies flip(bb_frequencies)] , [mn + se; flip(mn - se)]', colors(1,:), 'FaceAlpha', .5);
%     hold on;
%     plot(bb_frequencies, mn, 'o-', 'Color', colors(1,:), 'LineWidth', 2)
% 
%     
%     mn = lgmn(spectraBlank{dd},2);
%     se = std(spectraBlank{dd},[], 2) / sqrt(size(spectraBlank{dd},2));    
%     fill([bb_frequencies flip(bb_frequencies)] , [mn + se; flip(mn - se)]', colors(4,:), 'FaceAlpha', .5);    
%     plot(bb_frequencies, mn, '-o', 'Color', colors(4,:), 'LineWidth', 2)
%     
%     
%     set(gca, 'XScale', 'log', 'YScale', 'log', ...
%         'YLim', 10.^[1 2], 'XLim', [60 150], ...
%         'XTick', 60:12:150, 'XGrid', 'on', 'LineWidth', 2, 'FontSize', 20)
%     title(str{dd})
%     xlabel('Frequency (Hz)')
%     ylabel('Power (fT^2)')
% end
% 



fH = figure(dd); clf, set(gcf, 'Color', 'w')

xt = 60:12:160;
yl = 10.^[0.9 2.1];
for whichSubject = 1:8
    
    for dd = 1:2
        
        subplot(4,8, whichSubject + (dd-1)*8);
        cla;
        
        bootstat = bootstrp(100, @(x) lgmn(x,1), spectraFull{whichSubject, dd}');
        % mn = median(bootstat);
        se = prctile(bootstat, [16 84],1);
        
        fill([bb_frequencies flip(bb_frequencies)] , [se(1,:) flip(se(2,:))], colors(1,:), 'FaceAlpha', .5);
        hold on;
        %plot(bb_frequencies, mn, 'o-', 'Color', colors(1,:), 'LineWidth', 2)
        
        
        bootstat = bootstrp(100, @(x) lgmn(x,1), spectraBlank{whichSubject, dd}');
        % mn = median(bootstat);
        se = prctile(bootstat, [16 84],1);
        
        fill([bb_frequencies flip(bb_frequencies)] , [se(1,:) flip(se(2,:))], colors(4,:), 'FaceAlpha', .5);
        %plot(bb_frequencies, mn, 'o-', 'Color', colors(4,:), 'LineWidth', 2)
        
        
        set(gca, 'XScale', 'linear', 'YScale', 'log', ...
            'YLim', yl, 'XLim', [60 153], ...
            'XTick', xt, 'XGrid', 'on', 'LineWidth', 2, 'FontSize', 12)
        
        title(sprintf('S%d %s', whichSubject, str{dd}))
        xlabel('Frequency (Hz)')
        ylabel('Power (fT^2)')
        
       
    end
    
    for dd = 1:2
        
        subplot(4,8, whichSubject + (dd+1)*8);
        cla;
        
        bootfull  = bootstrp(100, @(x) lgmn(x,1), spectraFull{whichSubject, dd}');
        bootblank = bootstrp(100, @(x) lgmn(x,1), spectraBlank{whichSubject, dd}');
        
        bootdiff = (bootfull - bootblank) ./ bootblank;
        se = prctile(bootdiff, [16 84],1)*100;
        
        fill([bb_frequencies flip(bb_frequencies)] , [se(1,:) flip(se(2,:))], colors(1,:), 'FaceAlpha', .5);
        hold on, plot(bb_frequencies, 0*bb_frequencies, 'k--')
        
        set(gca, 'XScale', 'linear', 'YScale', 'linear', ...
            'YLim', 60*[-1 1], 'XLim', [60 150], ...
            'XTick', xt, 'XGrid', 'on', 'LineWidth', 2, 'FontSize', 12)
        
                title(str{dd})
        xlabel('Frequency (Hz)')
        ylabel('Percent above blank')

    end
    
end

 %% Outlier 
 
 figure, 
 for ii = 1:numSubjects
    subplot(3,3,ii) 
    histogram(spectraFull{ii,1}(11,:))
     
 end

 %% group average
 fH = figure; clf, set(gcf, 'Color', 'w')

 percentdiff = zeros(numSubjects, numFrequencies, 2);
 for dd = 1:2
     subplot(2,1,1), %cla;
     for whichSubject = 1:numSubjects                        
         fullmn  = lgmn(spectraFull{whichSubject, dd},2);
         blankmn = lgmn(spectraBlank{whichSubject, dd},2); 
         
         
         percentdiff(whichSubject,:,dd) = (fullmn - blankmn)./blankmn;
     end
     
              
     bootfull  = bootstrp(100, @(x) lgmn(x,1), spectraFull{whichSubject, dd}');
     bootblank = bootstrp(100, @(x) lgmn(x,1), spectraBlank{whichSubject, dd}');
     
     bootdiff = (bootfull - bootblank) ./ bootblank;
     se = prctile(bootdiff, [16 84],1)*100;
     
     fill([bb_frequencies flip(bb_frequencies)] , [se(1,:) flip(se(2,:))], colors(dd,:), 'FaceAlpha', .5);
     hold on, plot(bb_frequencies, 0*bb_frequencies, 'k--')
     
     set(gca, 'XScale', 'linear', 'YScale', 'linear', ...
         'YLim', 60*[-1 1], 'XLim', [60 150], ...
         'XTick', xt, 'XGrid', 'on', 'LineWidth', 2, 'FontSize', 12)
     
     title(str{dd})
     xlabel('Frequency (Hz)')
     ylabel('Percent above blank')

     
 end
   