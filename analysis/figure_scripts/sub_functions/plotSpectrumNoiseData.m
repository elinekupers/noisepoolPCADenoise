
figureDir       = fullfile(dfdRootPath, 'analysis', 'figures'); % Where to save images?
dataDir         = fullfile(dfdRootPath, 'analysis', 'data');    % Where to save data?


figure;
channelIdx = 1;
for whichSubject    = 1:8;     % Subject 99 has the synthetic dataset.

load(fullfile(dataDir,sprintf('S0%d_noisedata.mat',whichSubject)));
load(fullfile(dataDir,sprintf('S0%d_conditions.mat',whichSubject)));
load(sprintf(fullfile(dataDir, 's%02d_denoisedData_bb.mat'),whichSubject));

design = zeros(length(conditions), 3);
design(conditions==1,1) = 1; % condition 1 is full field
design(conditions==5,2) = 1; % condition 5 is left field
design(conditions==7,3) = 1; % condition 7 is right field

blanktrials = all(design==0,2);
blanktrials = blanktrials(~badEpochs);

fullfieldtrials = logical(design(~badEpochs,1));

% Define axes
f = (0:999);
xl = [2 150];


% compute spectrum
spec = abs(fft(squeeze(noisedata(channelIdx,:,:))))/size(noisedata,2)*2;

subplot(4,2,whichSubject);
% plot median
plot(f, mean(spec(:,fullfieldtrials),2),'b'); hold on;
plot(f, mean(spec(:,blanktrials),2),'k');

% format x and y axes
set(gca, 'XLim', xl, 'XScale', 'log', 'YScale','log');

% label figure, add stimulus harmonic lines, and make it look nice
xlabel('Frequency (Hz)'); ylabel('Power (fT^2)');

end


