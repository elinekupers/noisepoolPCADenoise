function [fH,ch] = megPlotMap(vals,clims,fH,cm,ttl)
% stylize plotOnEgi map the way we want it to look

if notDefined('cm'),     cm = 'jet';    end  % colormap
if notDefined('ttl');    ttl = '';      end  % title string
if notDefined('fH'),     fH = gcf;      end

fs = 14; % font size
renderer = 'zbuffer';

ssm_plotOnMesh(vals);
colormap(cm);
if ~notDefined('clims'), set(gca, 'CLim', clims); end
set(fH, 'Renderer', renderer);
tH = title(ttl);
set(tH, 'FontSize', fs);
ch = colorbar;

end

%% Subroutine
function ssm_plotOnMesh(sensor_data)

% clip data to 157 points
if length(sensor_data) > 157, sensor_data = sensor_data(1:157); end

cfg=[];
%cfg.interpolation = 'nearest';

cfg.layout = ft_prepare_layout(cfg, data_hdr);
cfg.style='straight';

%cfg.style='blank';
%cfg.electrodes ='numbers';
cfg.colorbar='yes';
cfg.maplimits='maxmin';
cfg.data = sensor_data';
cfg.data(isnan(cfg.data)) = nanmedian(sensor_data);

topoplot(cfg,cfg.data);

end



