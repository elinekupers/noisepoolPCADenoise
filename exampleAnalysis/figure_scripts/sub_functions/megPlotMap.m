function [fH,ch] = megPlotMap(vals,clims,fH,cm,ttl,data_hdr,cfg)

if notDefined('cm'),       cm = 'jet';    end  % colormap
if notDefined('ttl');      ttl = '';      end  % title string
if notDefined('fH'),       fH = gcf;      end
if notDefined('data_hdr'), data_hdr = []; end
if notDefined('cfg'),      cfg = []; end

fs = 14; % font size
renderer = 'zbuffer';

ssm_plotOnMesh(vals, data_hdr, cfg);
colormap(cm);
if ~notDefined('clims'), set(gca, 'CLim', clims); end
set(fH, 'Renderer', renderer);
tH = title(ttl);
set(tH, 'FontSize', fs);
ch = colorbar;

end

%% Subroutine
function ssm_plotOnMesh(sensor_data, data_hdr, cfg)

if notDefined('data_hdr'); data_hdr = load('meg160_example_hdr.mat');
    data_hdr = data_hdr.hdr;
else
    data_hdr = load(data_hdr); data_hdr = data_hdr.hdr;
end

if notDefined('cfg');
    % don't do anything. we don't need it
else
    cfg = load(cfg); cfg = cfg.cfg;
end

% Define plotting options
cfg.layout = ft_prepare_layout(cfg, data_hdr);
cfg.style='straight';
%cfg.interpolation = 'nearest';
%cfg.style='blank';
%cfg.electrodes ='numbers';
cfg.colorbar='yes';
cfg.maplimits='maxmin';

% clip data to 157 points
if length(sensor_data) < 157
    
    cfg.interpolation   = 'v4';
    cfg.data = sensor_data';
     topoplot(cfg,cfg.data);
%     ft_topoplotER(cfg,cfg.data);
    
elseif length(sensor_data) > 157,
    sensor_data = sensor_data(1:157);
    
    cfg.data = sensor_data';
    cfg.data(isnan(cfg.data)) = nanmedian(sensor_data);
    topoplot(cfg,cfg.data);
end

end



