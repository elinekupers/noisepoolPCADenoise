function freq = megGetSLandABfrequencies(f, T, slF)
% Get the stimulus-locked (sl_f), asynchronous broadband (ab_f), and keep
% (keep_f) frequencies and their indices (sl_i, ab_i, keep_i) into f
% (frequency vector). Broadband frequencies are all frequencies between min
% and max excluding harmonics of the stimulus-locked and line noise
% frequencies. Stimulus locked frequency is twice the stimulus flicker
% rate. Keep frequencies are the union of AB and SL freqencies.
%
% freq = megGetSLandABfrequencies(f, T, slF)
%
% Inputs:
%   f: vector of frequencies
%   T: epoch length (should be equal to 1/f(2))
%   slF: stimulus-locked frequency (twice the stimulus flicker rate)
% Output:
%   freq: structure containing the fields 'all', 'ab', 'sl', 'keep', 
%         'ab_i', 'sl_i', 'keep_i'

%% Check Inputs

% epoch length
if notDefined('T'), T = 1/f(2); end

% stimulus locked frequency
if notDefined('slF'), slF = 12/T; end

%% stimulus locked (2nd harmonic only).
% Should be same as slF, but this ensures that it is exactly a value in f
[~, sl_i] = min(abs(f - slF));
sl_f = f(sl_i);

%% Asyncrhonous Broadband - take all frequencies except: 
% (a) stimulus locked and harmonics,
% (b) line noise and harmonics
% (c) low frequencies (up to 7 Hz)

% drop the stimulus-locked frequncy and harmonics (frequencies within 3 Hz
% of multiples of the SL frequency, close to 15 Hz)
tmp     = (1:12)*slF;
%sl_drop  = sort([tmp-3 tmp-2 tmp-1 tmp tmp+1 tmp+2 tmp+3]);
sl_drop  = sort([tmp-1 tmp tmp+1]);

% drop line noise(frequencies within 3 Hz of multiples of 60 Hz)
tmp = (1:5) * 60;
%ln_drop   = sort([tmp-3 tmp-2 tmp-1 tmp tmp+1 tmp+2 tmp+3]);
ln_drop   = sort([tmp-1 tmp tmp+1]);

% low frequency drop
lf_drop = 0:60;

[~, ab_i]   = setdiff(round(f), [sl_drop ln_drop lf_drop]);
ab_f = f(ab_i);

%% keep frequencies - drop line noise and low frequencies, keep the rest
[~, keep_i]   = setdiff(round(f*T), [ln_drop lf_drop]);
keep_f = f(keep_i);

%% Return a struct containing frequencies and indices into these frequencies
freq.all    = f;        % all frequencies
freq.ab     = ab_f;     % frequencies used for computing broadband response
freq.sl     = sl_f;     % stimulus locked frequency (2*stimulus flicker rate)
freq.keep   = keep_f;   % all frequencies between min and max excluding line noise
freq.ab_i   = ab_i;     % index into f.all of ab frequencies
freq.sl_i   = sl_i;     % index into f.all of sl frequencies
freq.keep_i = keep_i;   % index into f.all of keep frequencies

return