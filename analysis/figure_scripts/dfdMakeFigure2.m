function dfdMakeFigure2()

%% Function to reproduce Figure 2 Methods figure 

% Eline Kupers, Helena X. Wang, Kaoru Amano, Kendrick N. Kay, David J.
% Heeger, Jonathan Winawer. (YEAR) Broadband spectral responses in visual
% cortex revealed by a new MEG denoising algorithm.
% (JOURNAL. VOLUME. ISSUE. DOI.)

whichSubject    = 1;     % Subject 1 has the example channel.
figureDir       = fullfile(dfdRootPath, 'analysis', 'figures'); % Where to save images?
dataDir         = fullfile(dfdRootPath, 'analysis', 'data');    % Where to save data?
saveFigures     = true;  % Save figures in the figure folder?
                                         
% Define plot colors
colors          = dfdGetColors(4);

[data,design,exampleIndex,exampleChannel] = prepareData(dataDir,whichSubject,4);

sensorData = data{1}; %take original data
design     = design{1}; %take corresponding design matrix

results = load(sprintf(fullfile(dataDir, 's%02d_denoisedData_bb.mat'),whichSubject));
badChannels = results.badChannels; clear results

%% Plot example time series - Fig 3A
% note that exampleIndex does not denote actual channel number, but number in
% the data matrix after removing bad channels (exampleChannel in actual
% channel space)

startEpoch = 85;%503;
nrEpochs   = 36;

chanNum0 = exampleIndex;
epochNum = startEpoch:startEpoch+nrEpochs-1;%1:36;
ts_time = (epochNum(1)+0.001):0.001:epochNum(end)+1;

fH = figure('position',[0,300,300,100*length(chanNum0)]);
for nch = 1:length(chanNum0)
    subplot(length(chanNum0),1,nch);
    ts = squeeze(sensorData(chanNum0(nch),:,epochNum));
    ts = ts(:);
    plot(ts_time,ts,'k');
    set(gca,'xtick',startEpoch:6:startEpoch+nrEpochs-1,'xlim',[startEpoch,startEpoch+nrEpochs-1],'ylim',3000*[-1,1],'xgrid','on');
    xlabel('Time (s)');
    makeprettyaxes(gca,9,9);
    % mark where 1 second of example epoch is
    %if nch==length(chanNum0), vline([3,4],'r'); end
end
% Print
if saveFigures
    figurewrite(fullfile(figureDir, 'figure2A_ts'),[],0,'.',1);
end

%% Plot power spectrum of the example time series - Fig 3B
% chanNum0 41, epoch 4

epochNum = startEpoch+4;
spec = abs(fft(squeeze(sensorData(chanNum0(end),:,:))))/size(sensorData,2)*2;
this_data = spec(:,epochNum).^2; % calculate power

f = (0:999);
xl = [8 150];
fok = f;
fok(f<=xl(1) | f>=xl(2) | mod(f,60) < 1 | mod(f,60) > 59) = [];

fH = figure('position',[0,300,300,200]); hold on;

plot(fok,this_data(fok+1,:),'k-');
xt = [12:12:144];
yt = 0:5; yl=[yt(1),yt(end)];
set(gca, 'XLim', [8 150], 'XTick', xt, 'XScale', 'log', ...
    'ytick',10.^yt, 'ylim',10.^yl,'YScale', 'log');
ss = 12; yl = get(gca, 'YLim');
for ii =ss:ss:180, plot([ii ii], yl, 'k--'); end
xlabel('Frequency (Hz)');
ylabel(sprintf('Power (%s)', 'fT^2'));
makeprettyaxes(gca,9,9);

% Print
if saveFigures
    figurewrite(fullfile(figureDir, 'figure2B_spec'),[],0,'.',1);
end

%% Plot broadband and stimulus locked simmary metric - Fig 3C
% chanNum0 41

% define frequencies and functions
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

% Define functions to define noise pool and signal of interest
funcs           = {@(x)getstimlocked(x,sl_freq_i), @(x)getbroadband(x,keep_frequencies,1000)};

% define x and y axes
epochNum = startEpoch:startEpoch+nrEpochs-1;
spects_time = [epochNum(1):epochNum(end)]-0.5;
% ylims = {[0,500],[0,80]};
ylims = {[0,500],[0,40]};

fH = figure('position',[0,300,300,200]);
for whichfun = 1:2
    subplot(2,1,whichfun);
    % get spectrum for this channel and these epochs, and apply function
    spects = funcs{whichfun}(sensorData(chanNum0(end),:,epochNum));
    % plot
    plot(spects_time,spects,'ko-','markersize',4);
    set(gca,'xtick',startEpoch:6:startEpoch+nrEpochs-1,'xlim',[startEpoch,startEpoch+nrEpochs-1],'xgrid','on','ylim',ylims{whichfun});
    makeprettyaxes(gca,9,9);
end

if saveFigures
    figurewrite(fullfile(figureDir, 'figure2C_sl'),[],0,'.',1);
end

%% Find channel location - Figure 3A, inset
% fprintf('original channel number is: %d\n', megGetOrigChannel(chanNum0,badChannels,false));
figure;
chanloc = ones(1,size(sensorData,1))*0.5;
chanloc(chanNum0(end))=1;
megPlotMap(to157chan(chanloc,~badChannels,'zeros'),[],[],'gray'); colorbar off

% for drawing the schematic as in Fig 3A inset, can also just make a
% uniform surface, and draw a circle around the channel number manually. To
% find where a particular channel is, be sure to set
% cfg.electrodes ='numbers'
% when calling topoplot. Also be sure to use the original channel number
% and not chanNum0

%% Plot design matrix - Figure 3D

% blank, full, right, left
condColors   = [0 0 0; dfdGetColors(3)];

fH = figure('position',[0,300,200,300]);
design_colored = design;
design_colored(:,2) = design_colored(:,2)*2;
design_colored(:,3) = design_colored(:,3)*3;

imagesc(design_colored); colormap(condColors)
makeprettyaxes(gca,9,9); axis off;

if saveFigures
    figurewrite(fullfile(figureDir, 'figure2D_design'),[],-1,'.',1);
end

%% calculate beta - Figure 3E
nepochs = size(sensorData,3);
nboot = 1000;
epochs_boot = randi(nepochs,nboot,nepochs);

% define figure properties
ylims = {[0,150],[0,5]};
fH = figure('position',[0,300,150,200]);

% mean subtract design matrix
design2 = bsxfun(@minus, design, mean(design));

for whichfun = 1:2 % loop through the functions (1: stimlocked, 2: broadband)
    beta = [];
    % calculate spectrum and substract mean
    spects = funcs{whichfun}(sensorData(chanNum0(end),:,:))';
    spects = bsxfun(@minus, spects, mean(spects));
    
    % each iteration we get 3 betas corresponding to the three conditions
    % do this for 100 bootstraps to obtain a distribution: [3 x 100] matrix
    for nn = 1:nboot
        curr_boot = epochs_boot(nn,:);
        curr_design = design2(curr_boot,:);
        curr_datast = spects(curr_boot,:);
        beta_boot = curr_design \ curr_datast;
        beta = cat(2,beta,beta_boot);
    end
    % calculate the median and std of bootstrapped beta
    beta_range = prctile(beta,[16 50 84],2)';
    beta_se = diff(beta_range([1 3],:))/2;
    
    % plot it and make it pretty
    subplot(2,1,whichfun);
    bar(beta_range(2,:),'EdgeColor','none','facecolor',[0.5,0.5,0.5]);
    errorbar2(1:3,beta_range(2,:),beta_range([1,3],:),1,'k-');
    set(gca,'xlim',[0.2,3.8],'ylim',ylims{whichfun},'ytick',ylims{whichfun});
    makeprettyaxes(gca,9,9);
end

if saveFigures
    figurewrite(fullfile(figureDir, 'figure2E_beta'),[],0,'.',1);
end