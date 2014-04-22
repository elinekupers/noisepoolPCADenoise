function [ab, sl, phases, bbPoly, slSNR] = megSummarizeSpectra(Amps, Ph, T, fmax)
% Convert a matrix of fourier coefficients (epochs x frequency) into
% vectors of EEG compoennts (stimulus-locked and broadband)
%
% [amps, phase, bbPoly] = eegSummarizeSpectra(Amps, Ph, T, f, fmax)
%
% Inputs:
%   Amps: Matrix of Fourier amplitudes, epoch x frequencies (x channel)
%   Ph:   Matrix of Fourier phases, epoch x frequencies (x channel)
%   T:    Epoch length, in seconds
%   fmax: maximum frequency used to compute broadband ECoG
%
% Outputs
%   amps: structure comprised of two vectors, ab (asyncrhonous broadband)
%           and sl (stimulus-locked), each (num epochs) x 1
%   phase: same as amps, but phase rather than amplitude. note that
%           phase.ab is an empty vector since the phase has no specific
%           meaning for the broadband signal
%   bbPoly:
%% CHECK INPUTS

% by default, use frequencies up to 150 Hz for broadband calculation
if notDefined('fmax'), fmax = 150; end

% if no Ph, make a matrix of NaN to avoid errors
if notDefined('Ph'), Ph = NaN(size(Amps)); end

% Get the stimulus-locked and broadband frequencies
freq = megGetSLandABfrequencies((0:fmax)/T, T, 12/T);

%% Asynchronous Broadband linear fit

% We assume a power law relationship between amplitude and frequency
% (amplitude is approximatley f^n, with n usually negative) as one component
% of the ECoG response. Because it is a power law, the relationship is
% linear in log-log space, meaning that log(amplitude) = n * log(frequency)
% + k. Because it is a power law, the linear relationship in log-log space
% holds whether the signal is defined as amplitude or squared amplitude
% (power). The only difference is that the exponent is doubled for power
% compared to amplitude. We fit the line in log-log space rather than
% fitting the power law to the untransformed spectrum because the variance
% in the log spectrum is approximately equal across frequencies.
%
% The fitting strategy is to take the power from many epochs, and fit a
% single linear function in log-log space to determine the slope (the power
% law exponent). We then find the intercept for individual epochs to get a
% measure of the broadband power from each epoch. Only a subset of
% frequencies are used for the fits; we exclude the frequencies near the
% even harmonics of the stimulus locked frequency (typicaly 15 Hz=2f) and
% near line noise and harmonics, as well as frequencies above fmax
% (typcially 150 Hz) or below fmin (typically 8 Hz).

ab = zeros(size(Amps,1), size(Amps,3));
bbPoly = zeros(2, size(Amps,3));
% loop over channels
for ii = 1:size(Amps,3)
    
    % fit the spectral data with a line in log-log space
    
    % Take the log of the frequencies and squared amplitude to get log-log space
    Y = log(Amps(:,freq.ab_i, ii).^2);
    X = repmat(log(freq.ab), size(Y, 1), 1);
    
    % fit ONE line to all the concatenated trials to get a single slope
    p = polyfit(X(:), Y(:), 1);
    
    % evaluate the fitted broadband line at the stimulus-locked frequency
    ssa = polyval(p, log(freq.sl));
    
    % calculate the residuals from the linear prediction to get AB reponses on
    % individual epochs
    residual = bsxfun(@minus, Y, polyval(p, log(freq.ab)));
    
    % take the mean of the residuals from each epoch and add the interecept
    % from the linear fit to get the  intercept for that epoch
    intercepts = nanmean(residual,2) + ssa;
    
    % broadband based on linear fit (we exponentiate to have units of power,
    % rather than log power)
    ab(:,ii) = exp(intercepts);
    
    % store the polynomial fits
    bbPoly(:,ii) = p;
    
end

%% Stimulus-locked time series (1st harmonic only)

sl     = squeeze(Amps(:,freq.sl_i,:));
phases = squeeze(Ph(:,freq.sl_i, :));

slnoise = squeeze(mean(Amps(:,[freq.sl_i-1 freq.sl_i+1],:),2));
slSNR  = sl ./ slnoise;

% other way of computing slSNR:
% sum [12,24,36,48] / sum [11,12,13,23,24,25,35,36,37,47,48,49]
% sl4h_i = (1:4)*(freq.sl_i-1)+1;    % make sure freq.all(sl4h_i) = [12,24,36,48]
% sl4hnoise_i = [sl4h_i-1, sl4h_i, sl4h_i+1];
% sl_sum1 = squeeze(sum(Amps(:,sl4h_i,:),2));
% sl_sum2 = squeeze(sum(Amps(:,sl4hnoise_i,:),2));
% slSNR  = sl_sum1 ./ sl_sum2;

% if sum(isnan(slSNR(:))) ~= 0
%     warning(sprintf('Nans found: electrode %d',find(isnan(slSNR(1,:)))));
%     slSNR(isnan(slSNR)) = eps;
%     ab(isnan(slSNR))    = eps;
% end

% % Subtract the AB from the Stimulus-Locked. We sqrt the AB because it was
% % computed as power, whereas SL is computed as amplitude
% sl = sl - sqrt(ab);


return