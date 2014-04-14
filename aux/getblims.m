function clims2 = getblims(evalout,prc)
% figure out the axis limits

clims2 = [];
for fh = 1:size(evalout,2)
    blims = cat(4,evalout(:,fh).beta);       % concat across pcs [n x channels x perms x pcs]
    clims = zeros(size(blims,1),2);          % n x [lower,upper]
    for whichbeta = 1:size(blims,1)
        blims2 = squeeze(blims(whichbeta,:,:,:)); % [channels x perms x pcs]
        blims2 = squeeze(mean(blims2,2));          % average across perms
        
        lower = prctile(min(blims2),prc(fh,1));     % min across channels for each pc, then prctile
        upper = prctile(max(blims2),prc(fh,2));     % max across channels for each pc, then prctile
        
        clims(whichbeta,:) = [-1, 1]*max(abs([lower, upper]));
    end
    clims2 = cat(3,clims2,clims); % n x [lower,upper] x nfuncs
end
