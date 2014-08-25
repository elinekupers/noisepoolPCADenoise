%% define and load example data set
sessionNum = 3;
[sessionName,megDataDir] = megGetDataPaths(sessionNum);
%inputDataDir = '/Volumes/HelenaBackup/denoisesuite/tmpmeg/';
%datafile = fullfile(inputDataDir,'inputdata',sprintf('%sb2',sessionName));
%disp(datafile);
%load(datafile,'sensorData','design','badChannels');
[sensorData, design, badChannels, conditionNames, okEpochs] ...
    = megLoadData(fullfile(megDataDir,sessionName),1:6);

figuredir = 'manuscript_figs/figure_glm';

%% Plot example time series - Fig 3A 
% note that chanNum0 does not denote actual channel number, but number in
% the data matrix after removing bad channels

%chanNum0  = [2,25,41];
chanNum0 = 41; 
epochNum = 1:36;
ts_time = 0.001:0.001:epochNum(end);

fH = figure('position',[0,300,300,100*length(chanNum0)]);
for nch = 1:length(chanNum0)
    subplot(length(chanNum0),1,nch);
    ts = squeeze(sensorData(chanNum0(nch),:,epochNum));
    ts = ts(:);
    plot(ts_time,ts,'k');
    set(gca,'xtick',0:6:36,'xlim',[0,36],'ylim',3000*[-1,1],'xgrid','on');
    xlabel('Time (s)'); 
    makeprettyaxes(gca,9,9);
    % mark where 1 second of example epoch is 
    if nch==length(chanNum0), vline([3,4],'r'); end
end
% Print
%figurewrite(fullfile(figuredir, 'figure_ts'),[],0,'.',1);

%% Plot power spectrum of the example time series - Fig 3B
% chanNum0 41, epoch 4

epochNum = 4;
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
%figurewrite(fullfile(figuredir, 'figure_spec'),[],0,'.',1);

%% Plot broadband and stimulus locked time series - Fig 3C
% chanNum0 41

% define functions 
T = 1; fmax = 150;
freq = megGetSLandABfrequencies((0:fmax)/T, T, 12/T);
funcs = {@(x)getstimlocked(x,freq), @(x)getbroadband(x,freq)};
% define x and y axes 
epochNum = 1:36;
spects_time = [1:epochNum(end)]-0.5;
ylims = {[0,500],[0,80]};

fH = figure('position',[0,300,300,200]);
for whichfun = 1:2
    subplot(2,1,whichfun);
    % get spectrum for this channel and these epochs, and apply function 
    spects = funcs{whichfun}(sensorData(chanNum0(end),:,epochNum));
    % plot 
    plot(spects_time,spects,'ko-','markersize',4);
    set(gca,'xtick',0:6:36,'xlim',[0,36],'xgrid','on','ylim',ylims{whichfun});
    makeprettyaxes(gca,9,9);
end

%figurewrite(fullfile(figuredir, 'figure_sl'),[],0,'.',1);

%% Find channel location - Figure 3A, inset
fprintf('original channel number is: %d\n', megGetOrigChannel(chanNum0,badChannels,false));

chanloc = zeros(1,size(sensorData,1));
chanloc(chanNum0(end))=1;
megPlotMap(to157chan(chanloc,~badChannels,'zeros')); colorbar off

% for drawing the schematic as in Fig 3A inset, can also just make a
% uniform surface, and draw a circle around the channel number manually. To
% find where a particular channel is, be sure to set 
% cfg.electrodes ='numbers'
% when calling topoplot. Also be sure to use the original channel number
% and not chanNum0

%% Plot design matrix - Figure 3D

% blank, full, right, left
condColors   = [0,0,0; 63, 121, 204; 228, 65, 69; 116,183,74]/255;

fH = figure('position',[0,300,200,300]);
design_colored = design;
design_colored(:,2) = design_colored(:,2)*2;
design_colored(:,3) = design_colored(:,3)*3;

imagesc(design_colored); colormap(condColors)
makeprettyaxes(gca,9,9); axis off;

%figurewrite(fullfile(figuredir, 'figure_design'),[],-1,'.',1);

%% calculate beta - Figure 3E
nepochs = size(sensorData,3);
nboot = 100;
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

%figurewrite(fullfile(figuredir, 'figure_beta'),[],0,'.',1);

%% find noise pool

% compute glm on stimulus-locked, using leave-one-out
design2 = bsxfun(@minus, design, mean(design));
datast = getstimlocked(sensorData,freq);
datast = bsxfun(@minus, datast, mean(datast));
epochs_test = (1:nepochs)';
% for each train/test permutation
modelfit = []; beta = []; r2perm = [];
for nn = 1:size(epochs_test,1)
    curr_test = epochs_test(nn,:);
    curr_train= setdiff(1:nepochs,curr_test);
    % do glm on training data
    % beta_train = [n x channels]
    beta_train= design2(curr_train,:)\ datast(curr_train,:);
    modelfit_test = design2(curr_test,:)*beta_train;
    % save prediction for this epcoh
    modelfit = cat(1,modelfit,modelfit_test); %[epochs x channels]
    beta     = cat(3,beta,    beta_train);    %[n x channels x perms]
end
% both data and predicted are [epochs x channels]; r2 = [1 x channels]
r2 = calccod(modelfit,datast(vectify(epochs_test'),:),[], 0);
[sortedr2,sortedinds] = sort(r2,'ascend');

%% Get and plot noise pool time series - Fig. 4B

nChan = 3; % number of channels to plot 
plot_hpf = true; % whether to high pass filter the data 

chanNum0  = sortedinds(1:nChan); % plot the ones with the worst fits
epochNum = 1:36;
ts_time = 0.001:0.001:epochNum(end);
if plot_hpf
    sensorDataFiltered = hpf(sensorData);
    data = sensorDataFiltered;
else
    data = sensorData;
end

fH = figure('position',[0,300,200,200]);
for nch = 1:length(chanNum0)
    subplot(length(chanNum0),1,nch);
    ts = squeeze(data(chanNum0(nch),:,epochNum));
    ts = ts(:);
    plot(ts_time,ts,'k');
    set(gca,'xtick',0:6:36,'xlim',[0,36],'ylim',500*[-1,1],'xgrid','on');
    axis off; %axis tight
    makeprettyaxes(gca,9,9);
    %if nch==length(chanNum0), vline(1:4,'r'); end
end

%figurewrite(fullfile(figuredir, 'figure_noisets'),[],0,'.',1);

%% Plot spatial map of noise pool on scalp - Fig. 4A
noisepool = false(1,size(sensorData,1));
noisepool(sortedinds(1:75))=true;
megPlotMap(to157chan(noisepool,~badChannels,'nans')); colorbar off

%figurewrite(fullfile(figuredir, 'figure_noisemap'),[],0,'.',1);

%% compute PCs and plot PC time series for 3 epochs  - Fig. 4C

for epochNum = 1:3
    % get noise pool
    noisedata = data(sortedinds(1:75),:,epochNum)';
    % make unit length
    temp = unitlengthfast(noisedata);
    % do svd
    [coef,u,eigvals] = princomp(temp);
    % make same size 
    pcs = bsxfun(@rdivide,u,std(u,[],1));
    
    % plot for each epoch  
    fH = figure('position',[0,100+epochNum*100,200,200]);
    for np = 1:3 % plot the first 3 PCs
        subplot(nChan,1,np);
        plot(1:1000,pcs(:,np),'k');
        set(gca,'ylim',4*[-1,1]);
        axis off;
        makeprettyaxes(gca,9,9);
    end
    suptitle(sprintf('Epoch %d',epochNum));
    %figurewrite(fullfile(figuredir,sprintf('figure_noisepcs_hpf_ep%d',epochNum)),[],0,'.',1);
end

%% plot the time series of three channnels - Fig. 4D
chanNum0  = sortedinds(end-nChan+1:end);
yoffset = 300;
colors = squeeze(varysat([0,0,1],[0.8,0.5,0.3]));

for epochNum = 1:3
    fH = figure('position',[0+epochNum*200,300,200,600]);
    for nch = 1:length(chanNum0)
        for npcs = [0,1,2]
            subplot(length(chanNum0),1,nch); hold on;
            ts = squeeze(data(chanNum0(nch),:,epochNum))';
            ts_denoised = ts - pcs(:,1:npcs)*(pcs(:,1:npcs)\ts);
            
            %plot(1:1000,ts,'k',1:1000,ts_denoised,'r');
            plot(1:1000,ts_denoised+200-yoffset*npcs,'color',colors(npcs+1,:));
            set(gca,'ylim',500*[-1,1]);
            axis off;
            makeprettyaxes(gca,9,9);
        end
        %figurewrite(fullfile(figuredir,sprintf('figure_denoisedts_hpf_ep%dpc%d',epochNum,npcs)),[],0,'.',0);
    end
    suptitle(sprintf('Epoch %d',epochNum));
end
