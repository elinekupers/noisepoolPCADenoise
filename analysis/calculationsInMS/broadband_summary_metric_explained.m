%% Simulate broadband in order to understand units

f = 1:200; % frequencies (Hz)
k = 10^8;  % Y intercept
n = -1;    % Slope in log-log space
response = k*f.^(n);

% add noise
response = response .* (.1+rand(size(f)));

bb_idx = setdiff(60:150, 12:12:200); % frequencies used to summarize broabband

% broadband summary metric 
getbb = @(x) geomean(x(bb_idx),2);
bb = getbb(response);

% linear fit in log log space
linpred = polyfit(log(f(bb_idx)), log(response(bb_idx)), 1);

% the crossover point is the temporal frequency at which the response
% matches the summary metric
crossover = geomean(bb_idx);
assert(isequal(geomean(bb_idx), crossover))

% visualize it
figure, plot(f, response,  ...
    f, exp(polyval(linpred, log(f))), 'r--', ...
    f([1 end]), [1 1]*bb)

set(gca, 'YScale', 'log', 'XScale', 'log');
yl = get(gca, 'YLim');
hold on
plot(crossover*[1 1], yl, 'k--')

legend('response', 'broadband summary metric', 'geometric mean of broadband frequencies')
