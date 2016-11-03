% Simulation of denoising on broadband signals

% Function handles
mnfrq   = @(x) mean(x,1);  % mean across frequencies (ie broadband power)
mntrial = @(x) mean(x,2);  % mean across trials 
se    = @(x,n) std(x,[],n)/sqrt(size(x,n)); % standard error of the mean



n  = 200;       % numnber of time points
nt = 100;       % number of trials
b  = 1;         % fixed component of broadband (same on all trials)
m  = 1;         % variable component of broadband (random across trials)
F  = NaN(n,nt); % power spectrum of undenoised data
FD = NaN(n,nt); % power spectrum of denoised data (variable component is reduced)
FB = NaN(n,nt); % power spectrum of baseline: data with smaller fixed component

% F: frequencies x trials

for ii = 1:nt
    trialnoise  = rand-.5;
    trialsignal = randn(1,n);

    s = trialsignal * (b+trialnoise*m);
    F(:,ii) = abs(fft(s));
    
    s = trialsignal * (b+trialnoise*m*.2);
    FD(:,ii) = abs(fft(s));
    
%     randwalk = randn(1,n);  
%     scalef = rand-.5;
 
    s = trialsignal * (b/2 + trialnoise*m);
    FB(:,ii) = abs(fft(s));
end

% Remove redundant fourier components
F = F(1:n/2, :);
FD = FD(1:n/2, :);
FB = FB(1:n/2, :);

% --- Plot spectra -----------------------------
f = 0:n/2-1;
figure(1), clf, set(gca, 'FontSize', 20, 'Color', .8 * [1 1 1]); hold on
plot(f, mntrial(F), 'b',  f, mntrial(FD), 'g', f, mntrial(FB), 'y', 'LineWidth', 2);

mn = mntrial(F)'; sem = se(F,2)';
fill([f flip(f)], [mn+sem flip(mn-sem)], 'b', 'FaceAlpha', 0.5)

mn = mntrial(FD)'; sem = se(FD,2)';
fill([f flip(f)], [mn+sem flip(mn-sem)], 'g', 'FaceAlpha', 0.5)

mn = mntrial(FB)'; sem = se(FB,2)';
fill([f flip(f)], [mn+sem flip(mn-sem)], 'k', 'FaceAlpha', 0.5)

plot(f, mntrial(F), 'b',  f, mntrial(FD), 'g', f, mntrial(FB), 'y', 'LineWidth', 2);

yl = get(gca, 'YLim');
set(gca, 'XScale', 'linear', 'YScale', 'linear', 'XLim', [0 n/2-1], 'YLim', [0 yl(2)]);
legend({'Undenoised', 'denoised', 'baseline'}, 'Location', 'Best')


% --- Plot broadband mean  -----------------------------

figure(2), clf, set(gca, 'FontSize', 20); hold on

mns   = mntrial([mnfrq(F); mnfrq(FD); mnfrq(FB)]);
seT   = se([mnfrq(F); mnfrq(FD); mnfrq(FB)],2);         % standard error over trials
seF   = se([mntrial(F) mntrial(FD) mntrial(FB)],1)';    % standard error over frequencies
h = bar([mns mns]'); hold on
hx = get(h, 'XData'); 
x = reshape([hx{:}], [], numel(h));
x = bsxfun(@plus, x, .225*[-1 0 1]);
errorbar(x,[mns mns]', [seT seF]', 'k.', 'LineWidth', 3)
set(gca, 'XTick', [1 2], 'XTickLabel', {'se(trials)', 'se(frq)'})
legend({'undenoised', 'denoised' 'baseline'}, 'Location', 'Best')
title('Broadband power: Mean±SEM')

% --- Plot broadband sem  -----------------------------

figure(3), clf, set(gca, 'FontSize', 20); hold on
seT   = se([mnfrq(F); mnfrq(FD); mnfrq(FB)],2);         % standard error over trials
seF   = se([mntrial(F) mntrial(FD) mntrial(FB)],1)';    % standard error over frequencies
bar([seT seF]'); hold on

set(gca, 'XTick', [1 2], 'XTickLabel', {'se across trials', 'se across frq'})
legend({'undenoised', 'denoised', 'baseline'}, 'Location', 'Best')
title('Broadband power: SEM')