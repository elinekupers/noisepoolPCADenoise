function fH = eegPlotMap(vals,fignum,cm,ttl,renderer,clims)
% stylize plotOnEgi map the way we want it to look

if notDefined('fignum'), fignum = []; end 
if notDefined('cm'),     cm = 'jet'; end  % colormap
if notDefined('ttl');    ttl = ''; end    % title string
if notDefined('renderer'), renderer = 'zbuffer'; end

fs = 18; % font size 

if ~isempty(fignum), fH = figure(fignum); clf; else fH = gcf; end
plotOnEgi(vals); 
colorbar; colormap(cm);
if ~notDefined('clims'), set(gca, 'CLim', clims); end
set(fH, 'Renderer', renderer);
tH = title(ttl);
set(tH, 'FontSize', fs)