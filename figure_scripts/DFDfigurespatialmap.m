
function [] = DFDfigurespatialmap(figNumbers, sessionNums,conditionNumbers,inputDataDir, fitDataStr, whichfun, figuredir,savefigures)
%% Reproduce requested spatial map figure
%  after using the 'Denoise Field Data' algorithm for the paper:
%
%   AUTHORS. YEAR. TITLE. JOURNAL. VOLUME. ISSUE. DOI.
%
% DFDfigurespatialmap(sessionNum,conditionNumbers,inputDataDir, ...
%                               fitDataStr, whichfun, figuredir,savefigures)
%
% Inputs:
%   figNumbers:         Which paper figures do you want to reproduce(5,8,9)
%                       [default=5];
%   sessionNums:        Vector of data sets (1 to 8)
%                       [default=1:8]
%   conditionNumbers:   Vector of conditions to use (1 to 6)
%                       [default=1:6]
%   inputDataDir:       String to define data directory
%                       [default=ullfile(DFDrootpath, 'data')]
%   fitDataStr:         String to define which input preproc data set to use
%                       [default='b2fr_hpf2_fit10p1k']
%   whichfun:           ???
%                       [default=1]
%   epochname:          String defining data directory to save data
%                       [default=fullfile(DFDrootpath, 'data')]
%   figuredir:          String to define where to save figures
%                       [fullfile(DFDrootpath,'figures')]
%   savefigures:        Boolean whether to save figures or not
%                       [default=false]
%
% Outputs:
%   One figure containing an example subject showing a stimulus signal,
%   broad before and after using the denoising algorithm, separate for the three
%   conditions. 
%
% Example: Denoise broadband of session 1 and save data
%   allResults = DFDDenoiseWrapper(1)
% 
% 



if notDefined('opt'),    opt = struct(); end
if ~isfield(opt,'sessionNums'),          opt.sessionNums = 1:8; end % Do all subjects
if ~isfield(opt,'conditionNumbers'),    opt.conditionNumbers = 1:6; end % Do all conditions
if ~isfield(opt,'inputDataDir'),        opt.inputDataDir = fullfile(DFDrootpath, 'data'); end % 
if ~isfield(opt,'fitDataStr'),          opt.fitDataStr = 'b2fr_hpf2_fit10p1k'; end % We need to automate this / link it to the opt. 
if ~isfield(opt,'whichfun'),            opt.whichfun = '1'; end
if ~isfield(opt,'figuredir'),           opt.figuredir = fullfile(DFDrootpath,'figures'); end
if ~isfield(opt,'savefigures'),         opt.savefigures = false; end


% Define example session 
sessionDir = DFDgetdatapaths(sessionNums,conditionNumbers,inputDataDir);
% stimulus locked glm results 
slresults = load(fullfile(inputDataDir,'savedProcData',sprintf('%sb2frSL_fit10p1k',sessionDir))); % Fix this hardcoded part!
slresults = slresults.results;
% broadband glm results 
thisfile = fullfile(inputDataDir,'savedProcData',sprintf('%s%s',sessionDir,fitDataStr));
disp(thisfile); load(thisfile,'results');
% load 'badChannels' matrix
datafile = fullfile(inputDataDir,'savedProcData',sprintf('%sb2',sessionDir));
disp(datafile); load(datafile,'badChannels');

%%
for ii = 1:length(figNumbers)
    switch figNumbers
        case 5 % Plot the subject the output head surface is in vector format. 
            figure('position',[1,600,1400,800]);
            for icond = 1:3
                % get stimulus-locked snr 
                sl_snr1 = getsignalnoise(slresults.origmodel(1),icond, 'SNR');
                %clims_sl = [0, max(sl_snr1)];
                clims_sl = [0,25.6723];
                % get broadband snr for before and after denoising 
                ab_snr1 = getsignalnoise(results.origmodel(whichfun),  icond, 'SNR');
                ab_snr2 = getsignalnoise(results.finalmodel(whichfun), icond, 'SNR');
                clims_ab = [0, max([ab_snr1, 12.4445])];
                %clims_ab = [0, max([ab_snr1, ab_snr2])];

                % convert back into 157-channel space 
                ab_snr1a = to157chan(ab_snr1,~badChannels,'nans');
                ab_snr2a = to157chan(ab_snr2,~badChannels,'nans');
                sl_snr1a = to157chan(sl_snr1,~badChannels,'nans');

                % plot spatial maps 
                subplot(3,3,(icond-1)*3+1)
                [~,ch] = megPlotMap(sl_snr1a,clims_sl,gcf,'jet','Stimulus Locked Original');
                makeprettyaxes(gca,9,9);
                makeprettyaxes(ch,9,9);

                subplot(3,3,(icond-1)*3+2)
                [~,ch] = megPlotMap(ab_snr1a,clims_ab,gcf,'jet','Original');
                makeprettyaxes(gca,9,9);
                makeprettyaxes(ch,9,9);

                subplot(3,3,(icond-1)*3+3)
                [~,ch] = megPlotMap(ab_snr2a,clims_ab,gcf,'jet',sprintf('Denoised PC %d',results.pcnum(whichfun)));
                makeprettyaxes(gca,9,9);
                makeprettyaxes(ch,9,9);
            end

            if savefigures
                figurewrite(fullfile(figuredir,'figure5_examplesubject3'),[],0,'.',1);
            end
            
            
        case 8 % Plot all subjects - full condition, post minus pre - Fig. 9
            % the output head surface is in vector format.
            figure('position',[1,600,1400,800]);
            icond = 1;
            for k = 1:length(sessionNums)
                sessionDir = DFDgetdatapaths(sessionNums(k),conditionNumbers, inputDataDir);
                thisfile = fullfile(inputDataDir,'savedProcData',sprintf('%s%s',sessionDir,fitDataStr));
                disp(thisfile); load(thisfile,'results');
                datafile = fullfile(inputDataDir,'savedProcData',sprintf('%sb2',sessionDir));
                disp(datafile); load(datafile,'badChannels');

                subplot(2,4,k);
                ab_snr1 = getsignalnoise(results.origmodel(whichfun),  icond, 'SNR');
                ab_snr2 = getsignalnoise(results.finalmodel(whichfun), icond, 'SNR');
                ab_snr_diff = to157chan(ab_snr2-ab_snr1,~badChannels,'nans');

                [~,ch] = megPlotMap(ab_snr_diff,[-10,10],gcf,jmaColors('coolhotcortex'));
                makeprettyaxes(gca,9,9);
                set(ch,'ytick',-10:5:10); 
                makeprettyaxes(ch,9,9); 
            end

            if savefigures
                figurewrite(fullfile(figuredir,'figure8_bbdiffall'),[],0,'.',1);
            end

        case 9  % All subjects - right minus left - Fig. 9b (before and after separately)
            figure('position',[1,600,1400,800]);
            %whichmodel = 'origmodel';
            whichmodel = 'finalmodel';
            for k = 1:length(sessionNums)
                sessionDir = DFDgetdatapaths(sessionNums(k),conditionNumbers, inputDataDir);
                thisfile = fullfile(inputDataDir,'savedProcData',sprintf('%s%s',sessionDir,fitDataStr));
                disp(thisfile); load(thisfile,'results');
                datafile = fullfile(inputDataDir,'savedProcData',sprintf('%sb2',sessionDir));
                disp(datafile); load(datafile,'badChannels');
                
                subplot(2,4,k);
            %     ab_snr1 = getsignalnoise(results.origmodel(whichfun), 2, 'SNR');
            %     ab_snr2 = getsignalnoise(results.finalmodel(whichfun),2, 'SNR');
            %     ab_snr_diff1 = to157chan(ab_snr2-ab_snr1,~badChannels,'nans');
            %     ab_snr1 = getsignalnoise(results.origmodel(whichfun), 3, 'SNR');
            %     ab_snr2 = getsignalnoise(results.finalmodel(whichfun),3, 'SNR');
            %     ab_snr_diff2 = to157chan(ab_snr2-ab_snr1,~badChannels,'nans');
            %     ab_snr_diff = ab_snr_diff1-ab_snr_diff2;
            
                ab_snr1 = getsignalnoise(results.(whichmodel)(whichfun), 2, 'SNR');
                ab_snr2 = getsignalnoise(results.(whichmodel)(whichfun), 3, 'SNR');
                ab_snr_diff = to157chan(ab_snr2-ab_snr1,~badChannels,'nans');
                
                [~,ch] = megPlotMap(ab_snr_diff,[-5,5],gcf,jmaColors('coolhotcortex'));
                makeprettyaxes(gca,9,9);
                %set(ch,'ytick',-20:10:20);
                set(ch,'ytick',-5:5:5);
                makeprettyaxes(ch,9,9); 
            end
            
            if savefigures
                figurewrite(fullfile(figuredir,'figure9ab_bbRightMLeft_after'),[],0,'.',1);
            end

    end
end

return