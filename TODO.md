1. Dependencies are:
	Statistics Toolbox (v 9.1)
	Signal Processing Toolbox (v 6.22)
	Neural Network Toolbox (v 8.2.1)
	
	Let's check whether there are many functions used. If just one or two, then maybe we can copy then into external or use replacement functions. 

2. hpf function has hardcoded numbers like sample rate is 1000. should we change this?

3. Save MEG data files as single sensor_data file rather than 6 files (one per condition)

4. Make figure scripts into functions rather than scripts

5. Figure 5: add 'full', 'left', 'right', and 'broadband' to titles

6. Figure out what whichFun means in DFDmakefigure6

7. Why is the script to make figure 6 saying that epochs are being grouped by 6?

8. Eliminate all suffices to denoised data files that are not necessary for clarification.
For example
	sensorDataStr      = 'b2';
	freqDefinition     = 'f';
	noisePoolSelection = 'r';
	
9. Eliminate conditionNumbers 

10. Deal with NaNs
	Option 1. Preprocessing:
	(a) Some epochs are labeled as bad for all channels. These should be eliminated from the dataset and the design matrix just before denoising (that is, delete rows from both arrays in preprocessing)
	(b) Some channels are labeled as bad for the entire experiment. Eliminate these before running the experiment (that is, eliminate columns from the data matr
	(c) For bad singletons (channels which have some but not all bad epochs), spatially interpolate across the array

	Option 2. Allow the denoisedata code to take optional inputs with bad channels and bad epochs
	(a) Allow opt.ignoreChannels
	(b) Allow opt.ignoreEpochs
	(c) For bad singletons, deal with in preprocessing

