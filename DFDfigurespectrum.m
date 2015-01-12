function DFDfigurespectrum(sessionNums, conditionNumbers, inputDataDir, fitDataStr, figuredir, savefigures)


%% define and load example data set

if notDefined('opt'),    opt = struct(); end
if ~isfield(opt,'sessionNums'),         opt.sessionNums = 1:8; end % Do all subjects
if ~isfield(opt,'conditionNumbers'),    opt.conditionNumbers = 1:6; end % Do all conditions
if ~isfield(opt,'inputDataDir'),        opt.inputDataDir = fullfile(DFDrootpath, 'data'); end % 
if ~isfield(opt,'fitDataStr'),          opt.fitDataStr = 'b2fr_hpf2_fit10p1k'; end % We need to automate this / link it to the opt. 
if ~isfield(opt,'figuredir'),           opt.figuredir = fullfile(DFDrootpath,'figures'); end
if ~isfield(opt,'savefigures'),         opt.savefigures = false; end


[sessionName,megDataDir] = DFDgetdatapaths(sessionNums,conditionNumbers, inputDataDir);
[sensorData, design, badChannels, ~, ~] ...
    = megLoadData(fullfile(megDataDir,sessionName),1:6);

%% Load fits and denoised data
% noisepool selection by SNR, highpass filtered, 75 channels in noisepool
% this file contains the denoisedts
datafile = fullfile(inputDataDir,'savedProcData',sprintf(['%s' opt.fitDataStr],sessionName));
disp(datafile); load(datafile);

%% Spectrum before and after
chanNum = 42;%26;
chanNum0 = megGetOrigChannel(chanNum,badChannels);

% %% check changes in beta and SNR
% results.origmodel(1).beta_md(:,chanNum0)
% results.finalmodel(1).beta_md(:,chanNum0)

% for k = 1:3
%     a = getsignalnoise(results.origmodel(1),k);
%     b = getsignalnoise(results.finalmodel(1),k);
%     [a(chanNum0),b(chanNum0)]
% end

%% plot spectrum - Figure 4AB

% define plot colors
colors = [63, 121, 204; 228, 65, 69; 116,183,74; 127,127,127]/255;

% full, right, left, off
condEpochs = {design(:,1)==1, design(:,2)==1, design(:,3)==1, all(design==0,2)};
% time domain data before and after denoising
data = {sensorData,denoisedts{1}};
% whether to average in log or not
avgLogFlg = false;

%% set up figure 4A
fH = figure('position',[0,300,500,500]); clf(fH);
% define axes
f = (0:999);
xl = [8 150];
fok = f;
fok(f<=xl(1) | f>=xl(2) | mod(f,60) < 2 | mod(f,60) > 58 ) = [];
xt = [12:12:72, 96,144];
yt = 1:5;
yl=[yt(1),yt(end)];

for dd = 1
    % compute spectrum
    spec = abs(fft(squeeze(data{dd}(chanNum0,:,:))))/size(data{dd},2)*2;

    hold on;
    for ii = [1,4]%1:length(condEpochs)
        % compute power for epochs corresponding to a condition and
        % trim data to specific frequencies
        this_data = spec(:,condEpochs{ii}).^2;
        this_data = this_data(fok+1,:);

        % compute median and confidence interval across epochs
        if avgLogFlg
            this_data_log = log10(this_data);
            %this_data_log(isinf(this_data_log)) = nan;
            mn = prctile(this_data_log,[16,50,84],2);
        else
            mn = prctile(this_data,[16,50,84],2);
        end
        mn(abs(mn-0)<0.001)= nan;

        % plot median
        plot(fok, mn(:,2),  '-',  'Color', colors(ii,:), 'LineWidth', 1);
        %plot(fok, mn(:,1),'Color', colors(ii,:));
        %plot(fok, mn(:,3),'Color', colors(ii,:));
    end

    % format x and y axes
    set(gca, 'XLim', xl, 'XTick', xt, 'XScale', 'log');
    if avgLogFlg
        set(gca,'ytick', yt, 'ylim',yl);
    else
        set(gca,'ytick',10.^yt, 'ylim',10.^yl,'YScale', 'log');
    end

    % label figure, add stimulus harmonic lines, and make it look nice
    xlabel('Frequency (Hz)'); ylabel('Power (fT^2)');
    title(sprintf('Channel %d', chanNum));
    yl2 = get(gca, 'YLim');
    for ii =12:12:180, plot([ii ii], yl2, 'k--'); end
    makeprettyaxes(gca,9,9);
end

if savefigures
    figurewrite(fullfile(figuredir,'figure4AFullSpectrumChannel42'),[],0,'.',1);
end





%% set up figure 4B
fH = figure('position',[0,300,200,350]); clf(fH);

% define axes
% f = (0:999);
%xl = [8 150];
xl = [60 150];
% fok = f;
% fok(f<=xl(1) | f>=xl(2) | mod(f,60) < 2 | mod(f,60) > 58 ) = [];
% xt = [12:12:72, 96,144];
yt = 1:2;
yl=[yt(1),yt(end)];

for dd = 1:2
    subplot(2,1,dd);
    
    % compute spectrum
    spec = abs(fft(squeeze(data{dd}(chanNum0,:,:))))/size(data{dd},2)*2;
    
    hold on;
    for ii = [1,4]%1:length(condEpochs)
        % compute power for epochs corresponding to a condition and
        % trim data to specific frequencies
        this_data = spec(:,condEpochs{ii}).^2;
        this_data = this_data(fok+1,:);
        
        % bootstrap for confidence intervals
        %         nboot = 1000;
        %         nepochs = size(this_data,2);
        %         epochs_boot = randi(nepochs,nboot,nepochs);
        %         epoch_vals = zeros(size(this_data,1),nboot);
        %         for nn = 1:nboot
        %             this_data_boot = this_data(:,epochs_boot(nn,:));
        %             if avgLogFlg
        %                 this_data_log = log10(this_data_boot);
        %                 epoch_vals(:,nn) = mean(this_data_log,2);
        %             else
        %                 epoch_vals(:,nn) = mean(this_data_boot,2);
        %             end
        %         end
        %         mn = prctile(epoch_vals,[2.5,50,97.5],2);
        
        % compute median and confidence interval across epochs
        if avgLogFlg
            this_data_log = log10(this_data);
            %this_data_log(isinf(this_data_log)) = nan;
            mn = prctile(this_data_log,[16,50,84],2);
        else
            mn = prctile(this_data,[16,50,84],2);
        end
        mn(abs(mn-0)<0.001)= nan;
        
        % plot median
        plot(fok, mn(:,2),  '-',  'Color', colors(ii,:), 'LineWidth', 1);
        %plot(fok, mn(:,1),'Color', colors(ii,:));
        %plot(fok, mn(:,3),'Color', colors(ii,:));
    end
    
    % format x and y axes
    set(gca, 'XLim', xl, 'XTick', xt, 'XScale', 'log');
    if avgLogFlg
        set(gca,'ytick', yt, 'ylim',yl);
    else
        set(gca,'ytick',10.^yt, 'ylim',10.^yl,'YScale', 'log');
    end
    
    % label figure, add stimulus harmonic lines, and make it look nice
    xlabel('Frequency (Hz)'); ylabel('Power (fT^2)');
    title(sprintf('Channel %d', chanNum));
    yl2 = get(gca, 'YLim');
    for ii =12:12:180, plot([ii ii], yl2, 'k--'); end
    makeprettyaxes(gca,9,9);
end

if savefigures
    figurewrite(fullfile(figuredir,'figure4bHighFreqSpectrumChannel42'),[],0,'.',1);
end

%% Bootstrap to get signal and noise - Fig 5C

f = (0:999);
fok = f;
fok(f<=60 | f>=150 ...
    | mod(f,60) < 2 | mod(f,60) > 58 ...
    | mod(f,12) < 2 | mod(f,12) > 10) = [];

meanXfreqs = {zeros(2,1000),zeros(2,1000)};
nboot = 1000;
for dd = 1:2 % for either pre and post denoising
    % compute spectrum
    spec = abs(fft(squeeze(data{dd}(chanNum0,:,:))))/size(data{dd},2)*2;
    
    % get power for full and blank epochs, at the specified frequencies
    this_data_full = spec(fok+1,condEpochs{1}).^2;
    this_data_blank = spec(fok+1,condEpochs{4}).^2;
    
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
end

% Set up figure and plot
fH = figure('position',[0,300,200,350]); clf(fH);
for dd = 1:2
    subplot(2,1,dd);
    %[n,x] = hist(meanXfreqs{dd}',30);
    %bar(x,n/1000,'barwidth',1,'edgecolor','none');
    %xlim([4,10]); set(gca,'xtick',4:2:10);
    vals = 22:0.2:42;
    n = hist(meanXfreqs{dd}',vals);
    bar(vals,n(:,1)/1000,'facecolor',colors(1,:),'edgecolor','none'); hold on;
    bar(vals,n(:,2)/1000,'facecolor',[0.5,0.5,0.5],'edgecolor','none');
    xlim([22,42]); set(gca,'xtick',22:4:42);
    ylim([0,0.45]);set(gca,'ytick',0:0.2:0.4);
    xlabel('Mean power (fT^2)'); ylabel('Fraction of bootstraps');
    makeprettyaxes(gca,9,9);
end

if savefigures
    figurewrite(fullfile(figuredir,'Figure4cFullDistributionBootDiff'),[],0,'.',1);
end