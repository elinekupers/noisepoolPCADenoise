%% Script to simulate the effect of denoising on the variance projected in

% One second time vector
t = 0:.001:1;

% Stimulus locked signal is sinusoid at 15 Hz
sl      = sin(2*pi*t*15);

% Noise is sampled from Gaussian distribution
noise   = randn(size(t));

% Data is a weighted sum of stimulus locked signal and noise
data    = sl*10+  .1* noise;

% Our principal component is also a weighted sum of the stimulus locked
% signal and noise. The weight on the stimulus locked signal is smaller in
% the pc than in the data. The noise has the same form as the noise in the
% data, but is a new sample.
noise   = randn(size(t));
pc      = sl+ .1*noise;

% Regress out the pc from the data
denoised = data - (data/pc) * pc;

% Plot
figure, set(gca, 'FontSize', 20); hold on
plot(t,data, t, pc, t, denoised)
legend({'data', 'pc', 'denoised data'}, 'Location', 'best')
xlim([0 .2])

% Note that the 'denoised' data has more 'noise' than the original data!
%   The noise we added is broadband. What happened is that when we
%   projected out the PC, the stimulus locked component matched between the
%   PC and the data, so that we effectively projected in noise from the PC.

% % We can quantify the broadband signal as the variance in the time series
% after regressing out the stimulus locked signal
bb = @(x) var(x - (x/sl) * sl);
disp([bb(data) bb(denoised)])


% Demonstrate that subtracting one noise component from another causes the
% variance to increase
a=rand(1,100);
b=rand(1,100);
var(a)
var(b)
var(a-b)