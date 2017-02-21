%% dfdMakeAllFigures
% 
% Eline Kupers, Helena X. Wang, Kaoru Amano, Kendrick N. Kay, David J.
% Heeger, Jonathan Winawer. (YEAR) A non-invasive, quantitative study of
% broadband spectral responses in human visual cortex
% (JOURNAL. VOLUME. ISSUE. DOI.)

% This is a script that will make all figures except 1, 4 (Methods) and 15
% (Discussion) from the MS: A non-invasive, quantitative study of
% broadband spectral responses in human visual cortex

% Eline Kupers, Feb 2017

% Fieldtrip toolbox is needed to plot topographic maps (see:
% http://www.fieldtriptoolbox.org/)
if isempty(which('ft_prepare_layout')), dfdAddFieldtripPath, end

dfdMakeFigure2; % Spectrum of channel 42, S1 before denoising
dfdMakeFigure3; % Topographic maps of SL and BB S1 and across subjects
dfdMakeFigure5; % Spectrum of channel 42, S1 before and after denoising
dfdMakeFigure6; % Effect of denoising on broadband SNR
dfdMakeFigure7; % Effect of denoising on broadband signal and noise
dfdMakeFigure8; % Topographic maps of BB before and BB after denoising, S1 and across subjects
dfdMakeFigure9; % Topographic maps of BB after denoising, both-hemifield stimulus, all individual subjects
dfdMakeFigure10; % Comparison of MEG Denoise and control analyses
dfdMakeFigure11; % Topographic maps and comparisons of MEG Denoise and other denosing algorithm
dfdMakeFigure12; % Effect of denoising on stimulus-locked SNR
dfdMakeFigure13; % Topographic maps and comparisons of Cinet data (MEG Denoise and TSSS algorithm)
dfdMakeFigure14; % Microsaccade histograms and topographic meshes with(out) microsaccades
dfdMakeFigure16; % Panels of preprocessing pipeline Methods figure


% Make supplementary figures
dfdMakeFigureSF1; % Topographic maps of SL and BB before denoising, all stimulus conditions, all individual subjects
dfdMakeFigureSF2; % Topographic maps of BB after denoising, all stimulus conditions, all individual subjects


