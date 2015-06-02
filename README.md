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

Data 		: Contains Folders with example data 
		sets and a folder with saved processed 
		data matrices, in order to reproduce 
		figures faster after denoising.

Experiments 	: Scripts and functions specific to 
		denoise MEG or EEG data with an on/off
		steady state stimulus paradigm.

External 	: General functions from other 
		toolboxes, repositories or researchers 

Figure scripts 	: Scripts to make figures 4-12 from  
		the manuscript.

Figures		: Folder where figures will be saved
		in .eps format, if requested.

Funcs		: Functions that specify how to define
		the ‘noise’ channel pool, signal of 
		interest, and filtering options.

Scripts		: Scripts to preload and shape the 
		data before denoising.

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

5) Denoise data by projecting out x PCs, compute output of interest (using evalfun), fit GLM on output response, cross validate

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
Recreate figure 4 from manuscript. In the Matlab prompt, type:

dfdAddPaths
dfdMakeFigure4()
