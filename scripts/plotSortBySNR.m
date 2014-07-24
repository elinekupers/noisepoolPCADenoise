function plotSortBySNR(results,whichfun,plotType,whichConds,sortByPre)
% function to plot denoising results sorted by SNR across channels
% results is struct from denoisedata
% rest of arguments are optional

if notDefined('whichfun'),   whichfun = 1;    end  % which function, if results contains more than one
if notDefined('plotType'),   plotType = 'snr'; end % snr, s, or n 
if notDefined('whichConds'), whichConds = 1:3; end % 1:3, or 1, 2, 3 (full, right, left)
if notDefined('sortByPre'),  sortByPre = true; end

% nconds x nchannels
ab_signal1 = abs(results.origmodel(whichfun).beta_md);
ab_noise1  = results.origmodel(whichfun).beta_se;
ab_signal2 = abs(results.finalmodel(whichfun).beta_md);
ab_noise2  = results.finalmodel(whichfun).beta_se;
ab_snr1    = ab_signal1./ab_noise1;
ab_snr2    = ab_signal2./ab_noise2;
nchan = size(ab_snr1,2);

% choose the values we want
switch plotType
    case 'snr'
        vals1 = max(ab_snr1(whichConds,:),[],1);
        vals2 = max(ab_snr2(whichConds,:),[],1);
    case 's'
        vals1 = max(ab_signal1(whichConds,:),[],1);
        vals2 = max(ab_signal2(whichConds,:),[],1);
    case 'n'
        vals1 = 1./mean(ab_noise1(whichConds,:),1);
        vals2 = 1./mean(ab_noise2(whichConds,:),1);
    case 'r2'
        vals1 = results.origmodel(whichfun).r2;
        vals2 = results.finalmodel(whichfun).r2;
end

% sort according to final snr 
[vals1_sorted,idx1] = sort(vals1);
[vals2_sorted,idx2] = sort(vals2);
np1 = results.noisepool(idx1);
np2 = results.noisepool(idx2);

if sortByPre  % colors go according to pre denoising 
    colors1 = jet(nchan);
    colors2 = colors1(idx2,:);
else          % colors go according to post denoising 
    colors2 = jet(nchan);
    colors1 = colors2(idx1,:);
end

figure('position',[200,100,700,800]), hold on;
for nchan = 1:nchan
    plot(nchan,vals1_sorted(nchan),'o','markersize',10,'color',getcolor(np1(nchan),colors1(nchan,:)));
    plot(nchan,vals2_sorted(nchan),'.','markersize',20,'color',getcolor(np2(nchan),colors2(nchan,:))); 
end
xlabel('channels');
ylabel(upper(plotType));
makeprettyaxes(gca,14,14); 


function color = getcolor(isnoise,c)
if isnoise
    color = [0.5,0.5,0.5];
else
    color = c;
end