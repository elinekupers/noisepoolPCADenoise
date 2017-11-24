
figureDir       = fullfile(dfdRootPath, 'analysis', 'figures'); % Where to save images?
dataDir         = fullfile(dfdRootPath, 'analysis', 'data');    % Where to save data?


figure;
% channelIdx = 1;
for whichSubject    = 1%:8;     % Subject 99 has the synthetic dataset.

    load(fullfile(dataDir,sprintf('S0%d_noisepcs.mat',whichSubject)));
% load(fullfile(dataDir,sprintf('S0%d_conditions.mat',whichSubject)));
% load(sprintf(fullfile(dataDir, 's%02d_denoisedData_bb.mat'),whichSubject));

figure; hold all;
for ii = 1:length(pcs)
    
    [~,~,~,~,explained] = pca(pcs{ii});
    plot(explained)
    
end

spec = fft(pcs{ii}()



end


