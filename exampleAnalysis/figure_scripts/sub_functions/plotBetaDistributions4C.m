function fH = plotBetaDistributions4C(data, exampleIndex, condEpochs, colors, saveFigures, figureDir)

%% Bootstrap to get signal and noise - Fig 5C

f = (0:999);
fok = f;
fok(f<=60 | f>=150 ...
    | mod(f,60) < 2 | mod(f,60) > 58 ...
    | mod(f,12) < 2 | mod(f,12) > 10) = [];

meanXfreqs = {zeros(2,1000),zeros(2,1000)};
nboot = 1000;
for exampleIndex = 1:157;
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
    end
    
    % Set up figure and plot
    figure(1); set(gcf,'position',[0,300,200,350]); clf;
    for dd = 1:2
        subplot(2,1,dd);
        %[n,x] = hist(meanXfreqs{dd}',30);
        %bar(x,n/1000,'barwidth',1,'edgecolor','none');
        %xlim([4,10]); set(gca,'xtick',4:2:10);
        vals = 18:0.2:50;
        n = hist(meanXfreqs{dd}',vals);
        bar(vals,n(:,1)/1000,'facecolor',colors(1,:),'edgecolor','none'); hold on;
        bar(vals,n(:,2)/1000,'facecolor',[0.5,0.5,0.5],'edgecolor','none');
        xlim([18,max(vals)]); set(gca,'xtick',18:4:max(vals));
        ylim([0,0.45]);set(gca,'ytick',0:0.2:0.4);
        xlabel('Mean power (fT^2)'); ylabel('Fraction of bootstraps');
        makeprettyaxes(gca,9,9);
    end
    
    if saveFigures
        figurewrite(fullfile(figureDir,sprintf('Figure4cFullDistributionBootDiff%d',exampleIndex)),[],0,'.',1);
    end
end
end