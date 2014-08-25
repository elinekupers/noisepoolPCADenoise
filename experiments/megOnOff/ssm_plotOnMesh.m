function fH = ssm_plotOnMesh(sensor_data, title_txt, figure_num, data_hdr, plotType)
%Plot a surface map of the MEG sensor data
% fH = ssm_plotOnMesh(sensor_data, [title_txt], [figure_num], meg_files, [plotType])
%
% Inputs
%   sensor_data: 1x157 vector of sensor data (if longer than 157, then
%           truncate)
%   title: optional string for title of plot
%   figure_num: optional integer to specify figure number
%   plotType: any of '2d', '3d', 'both'
%
% Output
%   fH: figure handle (integer)
%
% Example:
%   fH = ssm_plotOnMesh(randn(1,157), 'my random data', 999);


% check inputs
if ~exist('plotType', 'var') || isempty(plotType), plotType = 'both'; end

% check length of sensor data
if length(sensor_data) > 157, sensor_data = sensor_data(1:157); end



% Load the sensor positions in 3 space 
tmp = load('meg160xyz');
xyz = tmp.xyz;


%% Old mesh
% % Build an interpolant function of the z-dimension of the sensor locations
% F.z = TriScatteredInterp(xyz(:,1), xyz(:,2), xyz(:,3), 'natural');
% 
% % Build an interpolant function of the sensor data 
% F.data = TriScatteredInterp(xyz(:,1), xyz(:,2), sensor_data(1:157)',  'natural');
% 
% % Sample the interpolant at these grid points
% %   Grid ranges from -15 to +15 as the stored xyz values range from ~ -13 to 13
% ti = -15:1:15;
% [qx,qy] = meshgrid(ti,ti);
% qz = F.z(qx,qy);
% 
% sensor_data_interpolated = F.data(qx,qy);

%% New mesh - Requires fieldtrip

% add fieldtrip matlab code
if isempty(which('ft_analysispipeline')),
    addpath(genpath('/Volumes/server/Projects/MEG/code/fieldtrip'));
end

switch lower(plotType)
    case {'3d', 'both'}
        % set up figure
        if exist('figure_num', 'var'), fH = figure(figure_num); clf;
        else                           fH = figure; clf; end
        set(fH, 'Color', 'w')
        
        % Get labels
        colorbar
        ylabel('    Right       Left     ')
        xlabel('    Posterior       Anterior     ')
        zlabel('    Inferior       Superior     ')
        
        ft_plot_topo3d(xyz,sensor_data); hold on;
        label_add(xyz)
        
        % add a title if requested
        if exist('title_txt', 'var') && ~isempty(title_txt), title(title_txt); end
end

%% NEW TOPOPLOT

switch lower(plotType)
    case {'2d', 'both'}
        %rawdir = '/Volumes/server/Projects/MEG/SSMEG/02_SSMEG_02_28_2014/raw';
        %if notDefined('meg_files')
        %    meg_files = dir(fullfile(rawdir,'*.sqd'));
        %end
        %data_hdr = ft_read_header(fullfile(rawdir,meg_files(1).name));
        if notDefined('data_hdr')
            data_hdr = load('hdr'); data_hdr = data_hdr.hdr;
        end
        
        cfg=[];
        %cfg.interpolation = 'nearest';
        cfg.layout = ft_prepare_layout(cfg, data_hdr);
        cfg.style='straight';
        % cfg.style='blank';
        %cfg.electrodes ='numbers';
        cfg.colorbar='yes';
        cfg.maplimits='maxmin';
        cfg.data = sensor_data';

        
        for ii = 1:length(cfg.data)
            if isnan(cfg.data(ii))
                cfg.data(ii) = nanmedian(sensor_data);
            end
        end
        %figure; clf;
        topoplot(cfg,cfg.data);
        
        fH = gcf;
        
        % add a title if requested
        if exist('title_txt', 'var') && ~isempty(title_txt), title(title_txt); end
        
end

