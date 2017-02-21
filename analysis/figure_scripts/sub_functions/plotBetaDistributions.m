function fH = plotBetaDistributions(data, exampleIndex, exampleChannel, condEpochs, colors, saveFigures, figureDir, figNum, plotstr)

%% Bootstrap to get signal and noise - Fig 4C

f = (0:999);
fok = f;
fok(f<=60 | f>=150 ...
    | mod(f,60) < 2 | mod(f,60) > 58 ...
    | mod(f,12) < 2 | mod(f,12) > 10) = [];

meanXfreqs = {zeros(2,1000),zeros(2,1000)};
nboot = 1000;

for dd = 1:2 % for either pre and post denoising
    
    % compute spectrum
    spec = abs(fft(squeeze(data{dd}(exampleIndex,:,:))))/size(data{dd},2)*2;
    
    % get power for full and blank epochs, at the specified frequencies
    this_data_full = spec(fok+1,condEpochs{dd}{1}).^2;
    this_data_blank = spec(fok+1,condEpochs{dd}{4}).^2;
    
    % set up randomized indicies for bootstrapping
    nepochs_full = size(this_data_full,2);
    epochs_boot_full = randi(nepochs_full,nboot,nepochs_full);
    nepochs_blank = size(this_data_blank,2);
    epochs_boot_blank = randi(nepochs_blank,nboot,nepochs_blank);
    
    % bootstrap mean across epochs
    epoch_vals_full = zeros(length(fok),nboot);
    epoch_vals_blank = zeros(length(fok),nboot);
    for nn = 1:nboot
        this_data_full_boot  = this_data_full(:,epochs_boot_full(nn,:));
        this_data_blank_boot = this_data_blank(:,epochs_boot_blank(nn,:));
        epoch_vals_full(:,nn)  = mean(this_data_full_boot,2);
        epoch_vals_blank(:,nn) = mean(this_data_blank_boot,2);
    end
    
    % average across frequencies to get bootstrapped means for full or
    % blank conditions
    meanXfreqs{dd}(1,:) = mean(epoch_vals_full);
    meanXfreqs{dd}(2,:) = mean(epoch_vals_blank);
    
    pct = prctile(meanXfreqs{dd}(1,:), [16, 50, 84]);
    signal = pct(2);
    noise  = (pct(3)-pct(1))/2;
    fprintf('Signal: %4.2f\tNoise: %4.2f\tSNR:%4.2f\n', signal, noise, signal/noise);

    pct = prctile(meanXfreqs{dd}(2,:), [16, 50, 84]);
    signal = pct(2);
    noise  = (pct(3)-pct(1))/2;
    fprintf('Signal: %4.2f\tNoise: %4.2f\tSNR:%4.2f\n', signal, noise, signal/noise);

end

% Set up figure and plot
fH = figure; set(fH,'position',[0,200,400,400], 'Name', plotstr, 'NumberTitle', 'off');
if figNum == 2; dataToPlot = 1; elseif figNum == 5; dataToPlot = 1:2;
else disp('This figure number does not contain a distribution panel'); end
for dd = dataToPlot;
    subplot(2,2,dd);
    %[n,x] = hist(meanXfreqs{dd}',30);
    %bar(x,n/1000,'barwidth',1,'edgecolor','none');
    %xlim([4,10]); set(gca,'xtick',4:2:10);
    vals = 18:0.2:42;
    n = hist(meanXfreqs{dd}',vals);
    bar(vals,n(:,1)/1000,'facecolor',colors(1,:),'edgecolor','none'); hold on;
    bar(vals,n(:,2)/1000,'facecolor',[0.5,0.5,0.5],'edgecolor','none');
    xlim([min(vals),max(vals)]); set(gca,'xtick',min(vals):4:max(vals));
    ylim([0,0.5]);set(gca,'ytick',0:0.2:0.4);
    xlabel('Mean power (fT^2)'); ylabel('Fraction of bootstraps');
    makeprettyaxes(gca,9,9);
    
    % Make difference histogram
    subplot(2,2,dd+2);
    diffMeanFreq = meanXfreqs{dd}(1,:) - meanXfreqs{dd}(2,:);
    vals = 0:0.2:10;
    n = hist(diffMeanFreq',vals);
    bar(vals,n/1000,'facecolor','k','edgecolor','none'); hold on;
    xlim([0,10]); set(gca,'xtick',0:4:10);
    ylim([0,0.45]);set(gca,'ytick',0:0.2:0.4);
    
    xlabel('Mean power (fT^2)'); ylabel('Fraction of bootstraps');
    makeprettyaxes(gca,9,9);

    pct = prctile(diffMeanFreq, [16, 50, 84]);
    signal = pct(2);
    noise  = (pct(3)-pct(1))/2;
    fprintf('Signal: %4.2f\tNoise: %4.2f\tSNR:%4.2f\n', signal, noise, signal/noise);
    
    
    
end



if saveFigures
    if figNum == 2; fname = sprintf('Figure2cFullDistributionBootDiff%d',exampleChannel);
    elseif figNum == 5; fname = sprintf('Figure5bFullDistributionBootDiff%d',exampleChannel); end    
    figurewrite(fullfile(figureDir,fname),[],0,'.',1);    
end
end
