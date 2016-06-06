function fH = plotSpectraPanel4B(data, exampleIndex, condEpochs, avgLogFlg, colors, saveFigures, figureDir);

%% set up figure 4B
fH = figure; set(fH,'position',[0,300,200,350]);

% define axes
f = (0:999);
xl = [60 150];
fok = f;
fok(f<=xl(1) | f>=xl(2) | mod(f,60) < 2 | mod(f,60) > 58 ) = [];
xt = [12:12:72, 96,144];
yt = [1:2.5];
yl=[yt(1),yt(end)];

for exampleChannel = 1:157%42;
    clf;
    for dd = 1:2
        subplot(2,1,dd);
        
        % compute spectrum
        spec = abs(fft(squeeze(data{dd}(exampleChannel,:,:))))/size(data{dd},2)*2;
        
        hold on;
        for ii = [1,4]
            % compute power for epochs corresponding to a condition and
            % trim data to specific frequencies
            this_data = spec(:,condEpochs{dd}{ii}).^2;
            this_data = this_data(fok+1,:);
            
            % bootstrap for confidence intervals
            nboot = 1000;
            nepochs = size(this_data,2);
            epochs_boot = randi(nepochs,nboot,nepochs);
            epoch_vals = zeros(size(this_data,1),nboot);
            for nn = 1:nboot
                this_data_boot = this_data(:,epochs_boot(nn,:));
                if avgLogFlg
                    this_data_log = log10(this_data_boot);
                    epoch_vals(:,nn) = mean(this_data_log,2);
                else
                    epoch_vals(:,nn) = mean(this_data_boot,2);
                end
            end
            mn = prctile(epoch_vals,[2.5,50,97.5],2);
            
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
        title(sprintf('Channel %d', exampleChannel));
        yl2 = get(gca, 'YLim');
        for ii =12:12:180, plot([ii ii], yl2, 'k--'); end
        makeprettyaxes(gca,9,9);
    end
    
    if saveFigures
        figurewrite(sprintf(fullfile(figureDir,'figure4bHighFreqSpectrumChannel%d'),exampleChannel),[],0,'.',1);
    end
    
end