%% Script to reproduce Figure 4abc (Example data set) for one channel in example subject
%
% AUTHORS. TITLE. JOURNAL. YEAR.
%
% This figure will show (a) an the Fourier transformed amplitudes of an example 
% channel for baseline and full flicker condition. (b) High frequency
% portion before and after high pass filtering. (c) Distribution of broadband 
% and stimulus-locked SNR values.  
% each channel for the stimulus locked signal, broadband signals before and
% after using the denoising algorithm. The three separate conditions (Full,
% left, right hemifield stimulation are shown separately). 
%
% This script assumes that data is downloaded with the DFDdownloaddata
% function. 

%% Choices to make:

sessionNums         = 1;        % Choose 1:8 if you would like all the subjects.
conditionNumbers    = 1:6;      % Choose 1:6 to get all conditions: Full, 
                                % left, right (for 'on' 1,2,3 and 'off' 4,5,6)
                                % conditions. 
                                % If you like only a certain conditions e.g.
                                % full on and off define this variable as [1,4]
                                
inputDataDir = fullfile(DFDrootpath, 'data'); % Note: New data matrices will 
                                              % also get stored in the same folder.
                                              
fitDataStr          = 'b2fr_hpf2_fit10p1k'; % this string is needed to get the correct saved beta values
figureDir           = fullfile(DFDrootpath, 'figures');
saveFigures         = false;    % Save figures in the figure folder?

DFDfigurespectrum(sessionNums, conditionNumbers, inputDataDir, fitDataStr, figureDir, saveFigures)
