%% s_inspectNoiseData

% This script inspect the nature of the removed noise by plotting a bunch
% of figures for subject 1.

figureDir       = fullfile(dfdRootPath, 'analysis', 'figures'); % Where to save images?
dataDir         = fullfile(dfdRootPath, 'analysis', 'data');    % Where to save data?




for whichSubject    = 1%:8;
    
    load(fullfile(dataDir,sprintf('S0%d_noisedata.mat',whichSubject)));
    
    cmap = copper(size(noisedata,3));
    
    % Plot the variance explained by the PCs of the noise data
    
    explained = zeros(size(noisedata,1),size(noisedata,3));
    for ii = 1:size(noisedata,3)
        [coeff, score, latent, tsquared, explained(:,ii)] = pca(noisedata(:,:,ii)'); 
    end
    
    figure; hold on;
    plot(explained, 'Color',cmap(size(noisedata,3),:));
    plot(mean(explained,2),'k--','LineWidth',2)
    xlabel('Number of PC'); ylabel('Percentage of variance explained');
    title('Variance explained by PCs'); set(gca,'FontSize',20)
    
    
    % Plot spectrum of the noise data pcs
    load(fullfile(dataDir,sprintf('S0%d_noisepcs.mat',whichSubject)));
    
    figure; hold all;
    for ii = 1:size(pcs,1)
        amps(ii,:,:) = abs(fft(pcs{ii}(:,1:10))); 
%         plot(0:999,squeeze(amps(ii,:,:)));
    end
    
    for n_pcs = 1:10
    plot(1:1000, squeeze(mean(amps(:,:,n_pcs),1)),'k','LineWidth',1); end
    plot(1:1000, squeeze(mean(mean(amps,1),3)),'r','LineWidth',2);
    
    xlim([60 150]); ylim(10.^[1.5, 2.5]); set(gca,'XScale','log','YScale','log');
    xlabel('Frequency (Hz)'); ylabel('Amplitude (fT)');
    title(sprintf('Spectrum of noise 10 PCs'))
    
    
    % Plot mean spectrum of noise sensors
    figure; hold all;
    for jj = 1:75
        spec2 = fft(squeeze(noisedata(jj,:,:)));
        plot(0:999, mean(abs(spec2),2));
    end
    xlim([0 150]);   set(gca,'XScale','log','YScale','log');
    xlabel('Frequency (Hz)'); ylabel('Amplitude (fT)');
    title(sprintf('Mean spectrum of noisedata from %d noise sensors', 75))
    
    
    
    %% Inspect time series by plotting each pc time series separate
%     figure; 
%     for ii = 1:size(pcs)
%         clf;
%         for chan = 1:75
%             subplot(5,15,chan); title(sprintf('%d',chan));
%             plot(pcs{ii}(:,chan));
%         end
%      waitforbuttonpress;
%     end
    
    %% Or take one channel and compute the envelope
%     thischan = 5;
%     thisepoch = 10;
%     
%     thisnoise = noisedata(thischan,:,thisepoch);
%     envelope = hilbert(thisnoise);
%     figure; plot(0:999,thisnoise); hold on;
%     xlabel('Time (ms)')
%     ylabel('Amplitude (fT)')
%     title(sprintf('Time series of epoch %d, channel %d',thisepoch, thischan))
%     legend('Noise data', 'Envelope of noise data')
%     
%     figure; plot(0:999, abs(fft(envelope)));
%     xlim([0 150]); xlabel('Frequency (Hz)'); ylabel('Amplitude');
%     title('Fourier transform of one epoch time series envelope')
    
    
    %% Or the mean of envelopes across epochs
    
    envelope =  zeros(size(pcs{1},1),length(pcs));
    for ii = 1:867  
        envelope(:,ii) = hilbert(pcs{ii}(:,1));
    end
    
    amps_envelope = abs(fft(envelope))/length(size(envelope,1))*2;
    
    figure; plot(0:999,mean(amps_envelope,2));
    xlim([0 151]); ylim([200 450]); xlabel('Frequency (Hz)'); ylabel('Absolute amplitudes (fT)');
    title('Fourier transform of mean time series envelope');
    set(gca,'TickDir','out','FontSize',20); box off
    
end





