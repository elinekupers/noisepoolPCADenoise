function freq = megGetSLandABfrequencies(f, slF, lf_cutoff, tol)
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
%   slF: stimulus-locked frequency (twice the stimulus flicker rate)
%   lf_cutoff: exclude frequencies below this value
%   tol: exlude frequencies within tol of requested frequencies
%
% Output:
%   freq: structure containing the fields 'all', 'ab', 'sl', 'keep', 
%         'ab_i', 'sl_i', 'keep_i'
%
% Example:
%   freq = megGetSLandABfrequencies(0:199, 12, 60, 1.5);

%% Check Inputs

% stimulus locked frequency
if notDefined('slF'),       slF = [];       end
if notDefined('tol'),       tol = 1.5;      end
if notDefined('lf_cutoff'), lf_cutoff = 10; end

%% Asyncrhonous Broadband
% Exclude all frequencies that are close to a multiple of the
% stimulus-locked frequency
if isempty(slF), 
    sl_drop = []; 
else
    sl_drop  = f(mod(f, slF) <= tol | mod(f, slF) > slF - tol);
   
end

% Exclude all frequencies that are close to a multiple of the
% line noise frequency
ln_drop   = f(mod(f, 60) <= tol | mod(f, 60) > 60 - tol);

% Exclude all frequencies below lf_cutoff
lf_drop = f(f<lf_cutoff);

% Define the frequenies and indices into the frequencies used to compute
% broadband power
[ab_f, ab_i]   = setdiff(f, [sl_drop ln_drop lf_drop]);

%% Stimulus locked
if isempty(slF)
    sl_f = []; sl_i = [];
else
    [~, sl_i] = min(abs(f-slF)); sl_f = f(sl_i);
end

%% keep frequencies - for plotting
%   drop line noise and low frequencies, keep the rest (including those
%               used for broadband and stimulus locked compuation)
[keep_f, keep_i]   = setdiff(f, [ln_drop lf_drop]);

%% Return a struct containing frequencies and indices into these frequencies
freq.all    = f;        % all frequencies
freq.ab     = ab_f;     % frequencies used for computing broadband response
freq.sl     = sl_f;     % stimulus locked frequency 
freq.keep   = keep_f;   % all frequencies between min and max excluding line noise
freq.ab_i   = ab_i;     % index into f.all of ab frequencies
freq.sl_i   = sl_i;     % index into f.all of sl frequencies
freq.keep_i = keep_i;   % index into f.all of keep frequencies

return