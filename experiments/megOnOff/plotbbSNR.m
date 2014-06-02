function plotbbSNR(results,badChannels,whichbeta,whichfit,fH,type)

if notDefined('whichbeta'), whichbeta = 1:3; end
if notDefined('whichfit'),  whichfit  = 1;   end
if notDefined('fH')
    fH = figure('position',[1,600,1200 500]); set(fH, 'Color', 'w');
end
if notDefined('type'),     type = 'SNR'; end

switch type
    case 'S'
        ab_snr1 = max(abs(results.origmodel(whichfit).beta_md(whichbeta,:)),[],1);
        ab_snr2 = max(abs(results.finalmodel(whichfit).beta_md(whichbeta,:)),[],1);
    case 'N'
        ab_snr1 = mean(results.origmodel(whichfit).beta_se(whichbeta,:),1);
        ab_snr2 = mean(results.finalmodel(whichfit).beta_se(whichbeta,:),1);
    case 'SNR'  
        ab_snr1 = max(abs(results.origmodel(whichfit).beta_md(whichbeta,:)),[],1)./...
            mean(results.origmodel(whichfit).beta_se(whichbeta,:),1);
        ab_snr2 = max(abs(results.finalmodel(whichfit).beta_md(whichbeta,:)),[],1)./...
            mean(results.finalmodel(whichfit).beta_se(whichbeta,:),1);
end

clims_ab = [0, max([ab_snr1, ab_snr2])];
ab_snr1a = to157chan(ab_snr1,~badChannels,'nans');
ab_snr2a = to157chan(ab_snr2,~badChannels,'nans');

subplot(1,2,1);
megPlotMap(ab_snr1a,clims_ab,fH,'jet',sprintf('%s: Broadband Signal: Original',type));
subplot(1,2,2);
megPlotMap(ab_snr2a,clims_ab,fH,'jet',sprintf('%s: Broadband Signal: Denoised',type));