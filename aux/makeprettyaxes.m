function makeprettyaxes(ax,fs)
if ~exist('ax','var') || isempty(ax), ax = gca; end
if ~exist('fs','var') || isempty(fs), fs = 12; end
set(get(ax,'Title'),'FontSize',fs);
set(get(ax,'XLabel'),'FontSize',fs);
set(get(ax,'YLabel'),'FontSize',fs);
set(ax,'box','off','tickdir','out','ticklength',[0.025 0.025]);