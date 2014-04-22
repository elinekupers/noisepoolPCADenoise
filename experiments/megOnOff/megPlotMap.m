function fH = megPlotMap(vals,clims,fH,cm,ttl)
% stylize plotOnEgi map the way we want it to look

if notDefined('cm'),     cm = 'jet';    end  % colormap
if notDefined('ttl');    ttl = '';      end  % title string
if notDefined('fH'),     fH = gcf;      end

fs = 18; % font size 
renderer = 'zbuffer';

ssm_plotOnMesh(vals, [], [], [],'2d'); 
colorbar; colormap(cm);
if ~notDefined('clims'), set(gca, 'CLim', clims); end
set(fH, 'Renderer', renderer);
tH = title(ttl);
set(tH, 'FontSize', fs)