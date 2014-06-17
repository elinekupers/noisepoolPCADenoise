function plotbbSNR(results,badChannels,whichbeta,whichfit,fH,type)

if notDefined('whichbeta'), whichbeta = 1:3; end
if notDefined('whichfit'),  whichfit  = 1;   end
if notDefined('fH')
    fH = figure('position',[1,600,1200 500]); set(fH, 'Color', 'w');
end
if notDefined('type'),     type = 'SNR'; end

% get values 
ab_snr1 = getsignalnoise(results.origmodel(whichfit),whichbeta,type);
ab_snr2 = getsignalnoise(results.finalmodel(whichfit),whichbeta,type);

% map into original space 
ab_snr1a = to157chan(ab_snr1,~badChannels,'nans');
ab_snr2a = to157chan(ab_snr2,~badChannels,'nans');

% visualize
clims_ab = [0, max([ab_snr1, ab_snr2])];
subplot(1,2,1);
megPlotMap(ab_snr1a,clims_ab,fH,'jet',sprintf('%s: Original',type));
subplot(1,2,2);
megPlotMap(ab_snr2a,clims_ab,fH,'jet',sprintf('%s: Denoised, PC=%d',type,results.pcnum(whichfit)));