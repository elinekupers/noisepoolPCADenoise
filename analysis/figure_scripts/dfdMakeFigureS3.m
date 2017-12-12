function dfdMakeFigureS3()

%% Function to reproduce Supplementary Figure 3 (Summary statistics of removed noise timeseries)

% dfdMakeFigureS3()
%
% Eline Kupers, Helena X. Wang, Kaoru Amano, Kendrick N. Kay, David J.
% Heeger, Jonathan Winawer. (YEAR) A non-invasive, quantitative study of
% broadband spectral responses in human visual cortex
% (JOURNAL. VOLUME. ISSUE. DOI.)

% This function inspects the nature of the removed noise by plotting a bunch
% of figures with summary statistics for subject 1.

figureDir       = fullfile(dfdRootPath, 'analysis', 'figures'); % Where to save images?
dataDir         = fullfile(dfdRootPath, 'analysis', 'data');    % Where to save data?
saveFigures     = 1;


for whichSubject    = 1%:8;
    
    load(fullfile(dataDir,sprintf('S0%d_noisedata.mat',whichSubject)));
    
    cmap = copper(size(noisedata,3));
    
    % Plot the variance explained by the PCs of the noise data
    
    explained = zeros(size(noisedata,1),size(noisedata,3));
    
    for ii = 1:size(noisedata,3)
        [coeff, score, latent, tsquared, explained(:,ii)] = pca(noisedata(:,:,ii)');
    end
    
    
    figure; hold on;
    for ii = 1:size(noisedata,3); plot(explained(:,ii), 'Color', cmap(ii,:)); end
    plot(mean(explained,2),'k--','LineWidth',2)
    xlabel('Number of PC'); ylabel('Percentage of variance explained');
    title('Variance explained by PCs');
    makeprettyaxes(gca,9,9);
    
    if saveFigures
        figurewrite(fullfile(figureDir,'figureVarExpl'),[],0,'.',1);
    end
    
    figure; hold all;
    for ii = 1:size(noisedata,3); plot(explained(:,ii), 'Color', cmap(ii,:)); end
    plot(mean(explained,2),'k--','LineWidth',2);
    xlim([0 10]);
    xlabel('Number of PC'); ylabel('Percentage of variance explained');
    title('Variance explained by PCs');
    makeprettyaxes(gca,9,9);
    
    if saveFigures
        figurewrite(fullfile(figureDir,'figureVarExpl_inset'),[],0,'.',1);
    end
    
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
    
    xlim([60 150]); ylim(10.^[1.70, 2.1]); set(gca,'XScale','log','YScale','log');
    xlabel('Frequency (Hz)'); ylabel('Amplitude (fT)');
    title(sprintf('Spectrum of noise 10 PCs'));
    makeprettyaxes(gca,9,9);
    
    if saveFigures
        figurewrite(fullfile(figureDir,'figureSpectrumPCs'),[],0,'.',1);
    end
    
    
    %% Or take one channel and compute the envelope
    
    thispc = 1;
    thisepoch = 10;
    
    
    thisnoise = pcs{thisepoch}(:,thispc);
    envelope = hilbert(thisnoise);
    figure; plot(0:999,thisnoise); hold on;
    plot(0:999,abs(envelope));
    ylim([-4 4]);
    xlabel('Time (ms)')
    ylabel('Amplitude (fT)')
    title(sprintf('Time series of epoch %d, channel %d',thisepoch, thispc))
    makeprettyaxes(gca,9,9);
    
    legend('Noise data', 'Envelope of noise data')
    
    if saveFigures
        figurewrite(fullfile(figureDir,sprintf('figureNoiseData_pc%d_epoch%d',thispc,thisepoch)),[],0,'.',1);
    end
    
    figure; plot(0:999, abs(fft(abs(envelope))));
    xlim([0 150]);  ylim([0 1000]); xlabel('Frequency (Hz)'); ylabel('Amplitude');
    title('Fourier transform of one epoch time series envelope')
    makeprettyaxes(gca,9,9);
    
    if saveFigures
        figurewrite(fullfile(figureDir,sprintf('figureNoiseDataEnvelope_pc%d_epoch%d',thispc, thisepoch)),[],0,'.',1);
    end
    
    figure; plot(0:999, abs(fft(abs(envelope))));
%     ylim([0 275]);
    
    xlim([0 150]); xlabel('Frequency (Hz)'); ylabel('Amplitude');
    title('Fourier transform of one epoch time series envelope')
    makeprettyaxes(gca,9,9);
    
    if saveFigures
        figurewrite(fullfile(figureDir,sprintf('figureNoiseDataEnvelope_pc%d_epoch%d_inset',thispc, thisepoch)),[],0,'.',1);
    end
    
            
  
    
    
    
    %% Or the mean of envelopes across epochs
    
    envelope =  zeros(size(pcs{1},1),10,length(pcs));
    for ii = 1:867
        for jj = 1:10
            envelope(:,jj,ii) = hilbert(pcs{ii}(:,jj));
        end
    end
    
    envelope = reshape(envelope,1000,[]);
    amps_envelope = abs(fft(abs(envelope)));
    
    figure; plot(0:999,mean(amps_envelope,2), 'LineWidth',2);
    hold on; for ll = 1:7; plot(ll*[12 12], [0 250], 'k'); end
    xlim([0 151]); ylim([0 250]); xlabel('Frequency (Hz)'); ylabel('Absolute amplitudes (fT)');
    title('Fourier transform of mean time series envelope')
    makeprettyaxes(gca,9,9);
    
    if saveFigures
        figurewrite(fullfile(figureDir,'figureMeanNoiseDataEnvelope'),[],0,'.',1);
    end
    
end





