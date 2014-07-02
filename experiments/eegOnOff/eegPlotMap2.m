function fH = eegPlotMap2(vals,clims,fH,cm,ttl,renderer)
% stylize plotOnEgi map the way we want it to look

if notDefined('cm'),     cm = 'jet';    end  % colormap
if notDefined('ttl');    ttl = '';      end  % title string
if notDefined('fH'),     fH = gcf;      end
if notDefined('renderer'), renderer = 'zbuffer'; end

fs = 14; % font size 

plotOnEgi(vals); 
colorbar('location','eastoutside')
colormap(cm);
if ~notDefined('clims'), set(gca, 'CLim', clims); end
set(fH, 'Renderer', renderer);
tH = title(ttl);
set(tH, 'FontSize', fs);
