function pcchan = getTop10(results,whichfun)

if notDefined('whichfun'), whichfun = 1; end

% max across 3 conditions 
finalsnr = [getsignalnoise(results.origmodel(whichfun)); ...
    getsignalnoise(results.finalmodel(whichfun))];
% max across before and after 
finalsnr = max(finalsnr);
% exclude noise pool
finalsnr(results.noisepool) = -inf;
% sort 
[~,idx] = sort(finalsnr,'descend');
% find the top 10 
pcchan = false(size(results.noisepool));
pcchan(idx(1:10))= 1;