function fH = eegPlotLogSpectra(sensorData,condEpochs,chanNum,avgLogFlg,fH,lg)

if notDefined('avgLogFlg'), avgLogFlg = false; end
% make new figure, if no handle
if notDefined('fH'); fH = figure('Position',[0,600,700,500]); set(fH, 'Color', 'w'); end
if notDefined('lg'); lg = {'ON','OFF'}; end

spec = abs(fft(squeeze(sensorData(chanNum,:,:))))/size(sensorData,2)*2;

f = (0:999);
%   lower and upper bound of frequencies to plot (x lim)
xl = [8 150];
%   lower and upper bound of amplitudes to plot (y lim)
fok = f;
fok(f<=xl(1) | f>=xl(2) ...
    | mod(f,60) < 1 | mod(f,60) > 59 ...
    ... | mod(f,72) < 1 | mod(f,72) > 71 ...
    ... | abs(f-52) < 1 ...
    ) = [];

% plot colors
colors = [0.1 0.1 0.1; .6 .6 .6];
hold on;

for ii = 1:length(condEpochs)
    this_data = spec(:,condEpochs{ii}).^2;
    if avgLogFlg
        this_data_log = log10(this_data(fok+1,:));
        this_data_log(isinf(this_data_log)) = nan;
        plot(fok, nanmean(this_data_log,2), '-', 'Color', colors(ii,:), 'LineWidth', 2);
    else
        plot(fok, nanmean(this_data(fok+1,:),2),  '-',  'Color', colors(ii,:), 'LineWidth', 2);
    end
end

xt = [18:18:72,108,144];
yt = -2:3; yl=[yt(1),yt(end)];
set(gca, 'XLim', [8 150], 'XTick', xt, 'XScale', 'log', 'FontSize', 20);
if avgLogFlg
    set(gca,'ytick',yt, 'ylim',yl);
else
    set(gca,'ytick',10.^yt, 'ylim',10.^yl,'YScale', 'log');
end
xlabel('Frequency (Hz)');
ylabel(sprintf('Power (%s)', 'fT^2'));
title(sprintf('Channel %d', chanNum));
ss = 18; yl = get(gca, 'YLim');
for ii =ss:ss:180, plot([ii ii], yl, 'k--'); end
legend(lg);

