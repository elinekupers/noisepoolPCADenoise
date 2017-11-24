%% s_inspectNoiseData

% This script inspect the nature of the removed noise by plotting a bunch
% of figures for subject 1.

figureDir       = fullfile(dfdRootPath, 'analysis', 'figures'); % Where to save images?
dataDir         = fullfile(dfdRootPath, 'analysis', 'data');    % Where to save data?

nr_of_epochs_to_plot = 10; % or all epochs? =867
nr_of_channels_to_plot = 75;

cmap = copper(nr_of_epochs_to_plot);


for whichSubject    = 1%:8;
    
    load(fullfile(dataDir,sprintf('S0%d_noisedata.mat',whichSubject)));
    
    % Plot the variance explained by the PCs of the noise data
    figure; hold all;
    for ii = 1:nr_of_epochs_to_plot
        [~,~,~,~,explained] = pca(noisedata(:,:,ii));
        plot(explained, 'Color',cmap(ii,:));
    end
    xlabel('Number of PC'); ylabel('Percentage of variance explained');
    title(sprintf('Variance explained by PCs from the first %d epochs of %d noise sensors', nr_of_epochs_to_plot, nr_of_channels_to_plot))
    
    
    % Plot spectrum of the noise data pcs
    load(fullfile(dataDir,sprintf('S0%d_noisepcs.mat',whichSubject)));
    
    figure; hold all;
    for ii = 1:nr_of_epochs_to_plot
        spec = fft(pcs{ii});
        plot(0:999,abs(spec));
    end
    xlim([60 150]); ylim(10.^[0, 3]); set(gca,'XScale','log','YScale','log');
    xlabel('Frequency (Hz)'); ylabel('Amplitude (fT)');
    title(sprintf('Spectrum of PCs from the first %d epochs of %d noise sensors', nr_of_epochs_to_plot, nr_of_channels_to_plot))
    
    
    % Plot mean spectrum of noise sensors
    figure; hold all;
    for jj = 1:nr_of_channels_to_plot
        spec2 = fft(squeeze(noisedata(jj,:,:)));
        plot(0:999, mean(abs(spec2),2));
    end
    xlim([0 150]);   set(gca,'XScale','log','YScale','log');
    xlabel('Frequency (Hz)'); ylabel('Amplitude (fT)');
    title(sprintf('Mean spectrum of noisedata from %d noise sensors', nr_of_channels_to_plot))
    
    
    
    % Inspect time series
%     figure; 
%     for ii = 1:size(pcs)
%         clf;
%         for chan = 1:75
%             subplot(5,15,chan); title(sprintf('%d',chan));
%             plot(pcs{ii}(:,chan));
%         end
%      waitforbuttonpress;
%     end
    
    % Take one channel and compute the envelope
    
    thischan = 5;
    thisepoch = 10;
    
    thisnoise = noisedata(thischan,:,thisepoch);
    envelope = hilbert(thisnoise);
    figure; plot(0:999,thisnoise); hold on;
    plot(0:999,abs(envelope));
    
    figure; plot(0:999, abs(fft(envelope)));
    xlim([0 150]); xlabel('Frequency (Hz)')
    
    
end





