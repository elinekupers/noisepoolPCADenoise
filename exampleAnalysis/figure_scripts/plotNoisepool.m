whichSubject    = 9;        % Subject 1 is the example subject.
figureDir       = fullfile(dfdRootPath, 'exampleAnalysis', 'figures_rm1epoch_CiNet'); % Where to save images?
dataDir         = fullfile(dfdRootPath, 'exampleAnalysis', 'data');    % Where to save data?
saveFigures     = true;     % Save figures in the figure folder?
threshold       = 0;        % Set threshold for colormap. If no threshold set value to 0

if whichSubject > 8;
    data_hdr        = 'neuromag360_sample_hdr_combined.mat';
    cfg             = 'neuromag360_sample_cfg_combined.mat';
else
    data_hdr        = [];
end

%% Load denoised data of example subject
[data] = prepareData(dataDir,whichSubject,5);
bb = data{1};

% Get noisepool data
noisepool = double(bb.results.noisepool);

% Define color range
clims = [0,1];

% Plot it
[~,ch] = megPlotMap(noisepool,clims,gcf,'gray','Noisepool',data_hdr,cfg);
makeprettyaxes(gca,9,9);

% Save it
figurewrite(sprintf(fullfile(figureDir,'noisepool_subject%d_np40'),whichSubject),[],0,'.',1);
