readme.txt for saved processed data (saveProcData folder, in data folder)

Data set 4 has fitful75 with all 75 PCs used in the denoising

Other datasets with the fitull75 will have denoised results with datasets only up to 10 pcs (hack in denoisedata.m to save time). 
Yes, the name has to be changed, this will happen soon. 


From Helena's Readme file

Data: 
*m - data created by discarding beginning and end epochs of each stimulus presentations 
*b  - data created by allowing fewer 0 epochs, using a wider threshold (6 mm; used to be 3 mm) for selecting and averaging neighbors (8 mm for subj4/session5, subj5/session6, session7,session8, subj11/session10,subj12/session11)
*b2 - data created by allowing fewer 0 epochs, using the new data format (such that design matrix and condition order are as they occurred in the experiment)
*b3 - data created by allowing fewer 0 epochs, using the new data format, without physiological denoising

*hpf2 - high pass filtered with a sharp cutoff at 62 Hz and without stimulus harmonics
*p1k - bootstrapped 1000x rather than 100x

*epochGroup6  - denoising for 6 epochs at a time
*epochGroup6s- denoising for 6 epochs at a time, but shifted by 3 epochs, such that each denoising segment contains both ON and OFF data
*epochGroup6(s)o - makes sure that epoch groups are all the same lengths (discards those groups with fewer than 6 epochs because of bad epochs removed)

f  - uses new definition of ab_i for freq (includes more points; excludes 1 pt either side rather than 3)
r  - noise pool selection (and pc cutoff, if selected by algorithm) by SNR rather than by R2
lg - using getbroadbandlog as evalfun
fit10 - always use 10 PCs as cutoff