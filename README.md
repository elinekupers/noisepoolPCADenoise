General purpose denoising suite that denoises EEG/MEG/ECoG data

This denoising suite is dependent on MATLAB Version 8.4

——- With the following Matlab toolbox dependencies ---

Statistics Toolbox (v 9.1)
Signal Processing Toolbox (v 6.22)
Neural Network Toolbox (v 8.2.1)



---
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

