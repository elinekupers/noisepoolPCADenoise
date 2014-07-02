function newdata = hpfeeg(data)

% harmonics of of the stimulus locked
tmp = (1:14)*18; 
sl_drop  = sort([tmp-1 tmp tmp+1]);
% high pass filter with cutoff of 62 Hz, sharp cutoff, and excluding
% harmonics
newdata = filterdata(data,432,62,1,sl_drop);