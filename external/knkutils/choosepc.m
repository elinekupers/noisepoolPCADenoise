function chosen = choosepc(xvaltrend,pcchoose)

% this is the performance curve that starts at 0 (corresponding to 0 PCs)
curve = xvaltrend - xvaltrend(1);
% store the maximum of the curve
mx = max(curve);
% initialize (this will hold the best performance observed thus far)
best = -Inf;
for p=0:length(xvaltrend)-1
    % if better than best so far
    if curve(1+p) > best
        % record this number of PCs as the best
        chosen = p;
        best = curve(1+p);
        % if we are within opt.pcchoose of the max, then we stop.
        if best*pcchoose >= mx
        %if best >= mx * pcchoose
            break;
        end
    end
end
