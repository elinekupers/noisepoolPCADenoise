Welcome to our Denoise project code repository!

General purpose denoising suite that denoises EEG/MEG/ECoG data

This denoising suite is was developed on MATLAB Version 8.4

—————————————————————————————————————————————————————-
——-——-——-——-——- Matlab toolbox dependencies ——-———————
—————————————————————————————————————————————————————-

Statistics Toolbox (v 9.1)
Signal Processing Toolbox (v 6.22)
Neural Network Toolbox (v 8.2.1)

—————————————————————————————————————————————————————-
——-——-——-——-——-——- Folder structure ——-——-——-——-——-——-
—————————————————————————————————————————————————————-

alternativeAnalysis: (Under construction) Folder with 
		scripts to use alternative analyses, such as 
		ICA, to denoise MEG code.

exampleAnalysis:   	  Folder with code and possibility
		to download and store data to run analysis from 
		the paper and make all manuscript figures.
		
		Contains:
		- data:
		Empty folder to store data. Example data downloaded 
		by dfdDownloadSampleData will be written here by default.
		This folder contains a .gitignore file to prevent 
		large data files from being added to the repository. 
		- denoise_sample_data:
		Folder with specific code to denoise data of all subjects 
		for this particular steady state study.
		- figure_scripts:
		Functions to make figures 4-12 from the manuscript.
		- figures:
		 Folder where figures will be saved in .eps format,
		 if requested.
		- dfdDownloadSampleData.m: Function to download
		sample data for all subejcts from the web.

External 	: General functions from other 
		toolboxes, repositories or researchers 

denoisedata.m	: main function to denoise time series.

dfdDownloadData.m : Download example datasets from web.

dfdAddPaths.m	: Add paths with functions to make 
		this repository run smoothly.


—————————————————————————————————————————————————————-
——-——-——-——-——- General flow of scripts -——-——-——-——-—
—————————————————————————————————————————————————————-

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
——-——-——-——-——- Example -----------------——-——-——-——-—
—————————————————————————————————————————————————————-

% Prepare data sets.  In the Matlab prompt, type:
dfdAddPaths
dfdDownloadSampleData(1:8) 	% Slow. Do this once to download 8 data sets.
dfdDenoiseWrapper(1:8) 		% Slow. Do this once to denoise 8 sample data sets.

%  Recreate figure 4 from manuscript. 
dfdMakeFigure4()

%  Recreate all figures from manuscript. 
DFDmakeallfigures()

