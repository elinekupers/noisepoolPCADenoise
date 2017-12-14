%% dfdMakeAllFigures
% 
% Eline Kupers, Helena X. Wang, Kaoru Amano, Kendrick N. Kay, David J.
% Heeger, Jonathan Winawer. (YEAR) A non-invasive, quantitative study of
% broadband spectral responses in human visual cortex
% (JOURNAL. VOLUME. ISSUE. DOI.)

% This is a script that will make all figures except 1, 3 (Methods) and S7
% (Supplement) from the MS: A non-invasive, quantitative study of
% broadband spectral responses in human visual cortex

% Eline Kupers, Feb 2017

% Fieldtrip toolbox is needed to plot topographic maps (see:
% http://www.fieldtriptoolbox.org/)
if isempty(which('ft_prepare_layout')) 
    pth = '~/matlab/git/toolboxes/Fieldtrip/';
    dfdAddFieldtripPath(pth)
end

dfdMakeFigure2; % Panels of preprocessing pipeline Methods figure
dfdMakeFigure4; % Spectrum of channel 42, S1 before denoising
dfdMakeFigure5; % Topographic maps of SL and BB S1 and across subjects
dfdMakeFigure6; % Spectrum of channel 42, S1 before and after denoising
dfdMakeFigure7; % Effect of denoising on broadband SNR
dfdMakeFigure8; % Effect of denoising on broadband signal and noise
dfdMakeFigure9; % Topographic maps of BB before and BB after denoising, S1 and across subjects
dfdMakeFigure10; % Topographic maps of BB after denoising, both-hemifield stimulus, all individual subjects
dfdMakeFigure11; % Comparison of Noisepool-PCA and control analyses
dfdMakeFigure12; % Topographic maps and comparisons of Noisepool-PCA and other denosing algorithm
dfdMakeFigure13; % Topographic maps and comparisons of Cinet data (Noisepool-PCA and TSSS algorithm)
dfdMakeFigure14; % Microsaccade histograms and topographic meshes with(out) microsaccades

% Make supplementary figures
dfdMakeFigureS1; % Show the effect of denoising when varying number of channels in noisepool and pcs removed
dfdMakeFigureS2; % Show the effect of denoising when varying length of timepoints in epochs
dfdMakeFigureS3; % Synthetic dataset and the effect of denoising on known data
dfdMakeFigureS4; % Summary statistics of noise time series
dfdMakeFigureS5; % Topographic maps of SL and BB before denoising, all stimulus conditions, all individual subjects
dfdMakeFigureS6; % Topographic maps of BB after denoising, all stimulus conditions, all individual subjects
dfdMakeFigureS7; % Effect of denoising on stimulus-locked SNR


