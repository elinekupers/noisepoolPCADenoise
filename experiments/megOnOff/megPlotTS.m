function fH = megPlotTS(data,conds)

[fs,nepochs] = size(data);
ts     = (0:fs-1)/1000;
onsets = 0:nepochs-1;
colors = [.6 .6 .6; 0.1 0.1 0.1];

% Set up the figure 
fH = figure('Position',[0,600,800,500]);
set(fH, 'Color', 'w'); hold on
for ii = 1:nepochs
    plot(onsets(ii)+ts,data(:,ii),'Color', colors((conds(ii)==1)+1,:));
    if ii > 1 && conds(ii) ~= conds(ii-1)
        plot([onsets(ii),onsets(ii)],[-2000,2000],'k--');
    end
end

yl = std(data(:)) * [-5 5];
set(gca, 'XLim', [0,nepochs], 'XTick', 0:12:nepochs, ...
    'YLim', yl, 'YTick', -2000:500:2000, 'FontSize', 16);
xlabel('Time (s)')
ylabel('Strength (fT)')
