%% nppMakeAllFigures
% 
% Eline Kupers, Helena X. Wang, Kaoru Amano, Kendrick N. Kay, David J.
% Heeger, Jonathan Winawer. (2018) A non-invasive, quantitative study of
% broadband spectral responses in human visual cortex
% (PLOS ONE. VOLUME. ISSUE. DOI.)

% This is a script that will make all figures except 1 (Exp Design), 3
% (Methods) and S7 (Supplement) from the MS: A non-invasive, quantitative
% study of broadband spectral responses in human visual cortex

% Eline Kupers, Feb 2017

% Fieldtrip toolbox is needed to plot topographic maps (see:
% http://www.fieldtriptoolbox.org/)
if isempty(which('ft_prepare_layout')) 
    pth = '~/matlab/git/toolboxes/Fieldtrip/';
    nppAddFieldtripPath(pth)
end

nppMakeFigure2; % Panels of preprocessing pipeline Methods figure
nppMakeFigure4; % Spectrum of channel 42, S1 before denoising
nppMakeFigure5; % Topographic maps of SL and BB S1 and across subjects
nppMakeFigure6; % Spectrum of channel 42, S1 before and after denoising
nppMakeFigure7; % Effect of denoising on broadband SNR
nppMakeFigure8; % Effect of denoising on broadband signal and noise
nppMakeFigure9; % Topographic maps of BB before and BB after denoising, S1 and across subjects
nppMakeFigure10; % Topographic maps of BB after denoising, both-hemifield stimulus, all individual subjects
nppMakeFigure11; % Comparison of Noisepool-PCA and control analyses
nppMakeFigure12; % Topographic maps and comparisons of Noisepool-PCA and other denosing algorithm
nppMakeFigure13; % Topographic maps and comparisons of Cinet data (Noisepool-PCA and TSSS algorithm)
nppMakeFigure14; % Microsaccade histograms and topographic meshes with(out) microsaccades

% Make supplementary figures
nppMakeFigureS1; % Show the effect of denoising when varying number of channels in noisepool and pcs removed
nppMakeFigureS2; % Show the effect of denoising when varying length of timepoints in epochs
nppMakeFigureS3; % Synthetic dataset and the effect of denoising on known data
nppMakeFigureS4; % Summary statistics of noise time series
nppMakeFigureS5; % Topographic maps of SL and BB before denoising, all stimulus conditions, all individual subjects
nppMakeFigureS6; % Topographic maps of BB after denoising, all stimulus conditions, all individual subjects
nppMakeFigureS7; % Effect of denoising on stimulus-locked SNR


