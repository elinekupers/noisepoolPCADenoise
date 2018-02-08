Welcome to our Noisepool-PCA code repository!

General purpose denoising suite that denoises EEG/MEG/ECoG data

This denoising suite was developed on MATLAB Version 8.4 and is described in the manuscript in preparation,
*A non-invasive, quantitative study of broadband spectral responses in human visual cortex*
Eline Kupers, Helena X. Wang, Kaoru Amano, Kendrick N. Kay, David J. Heeger, Jonathan Winawer


**Matlab toolbox dependencies**
* Statistics Toolbox (v 9.1)
* Signal Processing Toolbox (v 6.22)
* Neural Network Toolbox (v 8.2.1)

**Other toolbox dependencies**
* Fieldtrip toolbox (v 20170102) - For visualization
* MGL (v 2.0) - For eyetracking analysis
* MrTools (v unknown) - For eyetracking analysis

**Folder structure**
* denoisedata.m:  Main function to denoise time series.
* nppAddPaths.m:  Add paths with functions to make this repository run smoothly.		
* nppAddFieldtripPaths.m: Add paths of Fieldtrip toolbox to plot data on a sensormap or to use ICA function.
* external (folder): Contains functions from other toolboxes, repositories or researchers 
* analysis (folder): Contains code to download, store, and analyze data from manuscript	
 * data (folder): Empty folder to store data. Example data downloaded by nppDownloadSampleData will be written here by default. This folder contains a .gitignore file to prevent large data files from being added to the repository. 
 * denoise_subfunctions (folder): Folder with subfunctions used by the nppDenoiseWrapper in order to denoise the data of all subjects for this particular steady state study.
 * figure_scripts (folder): Functions to make figures 4-14 for the manuscript.
 * figures (folder): Empty folder where figures made by figure_scipts will be saved. This folder contains a .gitignore file to prevent large image files from being added to the repository.
 * nppDownloadSampleData.m: Function to download sample data for all subjects from the web.
 * nppDenoiseWrapper.m: Function to denoise sample data for all subjects 
 * nppMakeAllFigures.m: Script to make all figures from the manuscript.

**General flow of denoise data function**


INPUT:

1) Data (channel x time x epoch)

2) Design matrix

3) Function to compute evoked response (evokedfun)

4) Function(s) of interest (evalfun)

---
WHAT THE MAIN FUNCTION COMPUTES:

1) Compute evoked response of input data (using evokedfun and data)

2) Fit GLM on evoked response and cross validate for each channel (using design matrix)

3) Select n channels as noise pool based on some criterion (e.g., R^2 of fits)

4) Compute PCA on noise pool

5) Denoise data by projecting out x PCs, compute output of interest (using evalfun), 
	fit GLM on output response, cross validate

6) Repeat for x+1 PCs, until we've tried some reasonable number of PCs

7) Select optimal number of PCs

---
OUTPUT:

1) Final model (GLM solution denoised with the optimal number of PCs)

2) All Fits with different number of PCs (mostly for posthoc analyses)

3) Noise pool (a vector of booleans)

4) Denoised output data (output of evalfun)

—————————————————————————————————————————————————————-
——-—- Example 1: Download raw data and denoise -------
—————————————————————————————————————————————————————-

	% Prepare data sets.  In the Matlab prompt, type:
	addpath(genpath(pwd))
	nppAddFieldtripPath

	% Download the sample data. Slow. Do this once to download 8 raw data sets.
	nppDownloadSampleData([],1:8,'raw');  

	% Denoise with exactly 10 PCs. Slow. Do this once to denoise 8 data sets.
	nppDenoiseWrapper(1:8,1); % 

	% Denoise up to 10 PCs. Slow. Do this once to denoise 8 data sets.
	nppDenoiseWrapper(1:8,2); 

	% Denoise with control methods. Slow. Do this once to denoise 8 data sets.
	nppDenoiseWrapper(1:8,3) 				 % 	

	%  Recreate figure 4 from manuscript. 
	nppMakeFigure4();

	%  Recreate all figures from manuscript. 
	nppMakeallfigures();

—————————————————————————————————————————————————————-
—— Example 2: Download denoised data and plot fig7 ---
—————————————————————————————————————————————————————-

	% Prepare data sets.  
	addpath(genpath(pwd));
	nppAddFieldtripPath();

	% Download sample data. Slow. Do this once to download 8 raw data sets.
	nppDownloadSampleData([],1:8,'denoised 10 pcs') 
	
	% Denoise with exactly 10 PCs. Slow. Do this once to denoise 8 data sets.
	nppDenoiseWrapper(1:8,1) 				

	%  Recreate figure 7 from manuscript. 
	nppMakeFigure7()
