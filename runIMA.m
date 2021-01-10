% performs time-frequency decomposition, PCA and runs Independent Modulator
% Analysis (IMA) on the matrix PC x spectra*ICs 
%
%
%[IMAweights,IMAsphere, meanpwr, freqvec, timevec,...
%    pcs, Trials, eigvec, pc, pwr, meanpwrCond,...
%    CondTimewindows, CondTrials, CondTimevec] = runIMA(eegdata,frames, tlimits, srate, varargin);
%
%
%
% Author: Johanna Wagner, Swartz Center for Computational Neuroscience, UC San Diego, 2019
% adapted from a function written by Julie Onton
%
%
%
% Example
% >>  [IMAweights,IMAsphere, meanpwr, freqvec, n_trials, ntw_trials,...
%    pcs, eigvec, pc, timefreq,...
%    meanpwrCond] = runIMA(EEG.icaact,EEG.pnts, [EEG.xmin EEG.xmax]*1000, EEG.srate,...
%    'plotcomps', [1 3 4 5 7 8 9 12 15 22], 'cycles', [6 0.5], 'timesout', 400,...
%     'frqlim', [3 120], 'freqscale', ''log', 'winsize', 1000, 'icatype', 'amica');   
%
%
%
% INPUTS:
% eegdata -- EEG data in 3D - components x time x trials can be either a single EEG data matrix or a cell array of two or more EEG data matrices
% plotcomps -- components to run IMA on
% frqlim -- [minfrq maxfrq] minimum and maximum frequencies to include in spectral decomposition
% freqscale -- ['linear',or 'log']
% pcfac -- [integer] regulates the dimension that the time windows of the spectral data will
%                      be reduced to before IMA - the smaller pcfac, the more dimensions will be
%                      retained ndims = sqrt((freqs*ICs)/pcfac) where freqs is the number of estimated
%                      frequencies and ICs is the number of ICs (default is
%                      7). This is to regulate the number of rows in
%                      proportion to the number of columns that go into ICA
% icatype -- ['string']  which ICA algorithm to run  {'amica' 'infomax'}  (default  'infomax')   
%
%% time frequency decomposition parameters see also newtimef()
% cycles -- [real] indicates the number of cycles for the time-frequency
%                          decomposition {default: 0}
%                       If 0, use FFTs and Hanning window tapering.
%                       If [real positive scalar], the number of cycles in each Morlet
%                          wavelet, held constant across frequencies.
%                       If [cycles cycles(2)] wavelet cycles increase with
%                          frequency beginning at cycles(1) and, if cycles(2) > 1,
%                          increasing to cycles(2) at the upper frequency,
%                        If cycles(2) = 0, use same window size for all frequencies
%                          (similar to FFT when cycles(1) = 1)
%                        If cycles(2) = 1, cycles do not increase (same as giving
%                           only one value for 'cycles'). This corresponds to a pure
%                           wavelet decomposition, same number of cycles at each frequency.
%                        If 0 < cycles(2) < 1, cycles increase linearly with frequency:
%                           from 0 --> FFT (same window width at all frequencies)
%                           to 1 --> wavelet (same number of cycles at all frequencies).
%                       The exact number of cycles in the highest frequency window is
%                       indicated in the command line output. Typical value: 'cycles', [3 0.5]
% timesout -- Number of output times (int<frames-winframes). Enter a
%                       negative value [-S] to subsample original times by S.
%                       Enter an array to obtain spectral decomposition at
%                       specific times (Note: The algorithm finds the closest time
%                       point in data; this could give a slightly unevenly spaced
%                       time array
% nfreqs -- number of output frequencies. For FFT, closest computed
%                       frequency will be returned. Overwrite 'padratio' effects
%                       for wavelets. {default: use 'padratio'}
% padratio --           FFT-length/winframes (2^k)                    {default: 4}
%                       Multiplies the number of output frequencies by dividing
%                       their spacing (standard FFT padding). When cycles~=0,
%                       frequency spacing is divided by padratio.
% wletmethod --['dftfilt'|'dftfilt2'|'dftfilt3'] Wavelet type to use.
%                       'dftfilt2' -> Morlet-variant wavelets, or Hanning DFT.
%                       'dftfilt3' -> Morlet wavelets.  See the timefreq() function
%                        for more details {default: 'dftfilt3'}
% winsize  --   [integer] in ms size of timewindows to use for time-frequency decomposition
%                         needs to be smaller than epochsize. Note: this parameter 
%                        is overwritten when the minimum frequency and number of cycles requires
%                        a longer time window. {default: cycles*(length of a cycle at lowest freq in sec)}
%
%
%
% OUTPUTS
% IMAweights - weights of IMA decomposition
% IMAsphere - spheres of IMA decomposition
% meanpwr - mean power spectra of single ICs
% freqvec - frequency vector
% pcs - number of dimensions used for PCA before IMA
% eigvec - pc backprojection in time
% pc - pc spectral backprojection
% pwr - spectral timefrequency decomposition (spectograms for each IC)
% meanpwrCond - mean power for each condition if run on STUDY otherwise
%              empty
% n_trials - number of trials in each condition
% ntw_trials - number of timewindows per trial



function [IMAweights,IMAsphere, meanpwr, freqvec ,n_trials, ntw_trial,...
    pcs, eigvec, pc, pwr, meanpwrCond] = runIMA(eegdata,frames, tlimits, srate, varargin);

if length(eegdata) > 1,
    tempdat = eegdata{1};
else
    tempdat = eegdata;
end

g = finputcheck(varargin, { 'plotcomps'     'integer'   []             [1:size(tempdat,1)]; ...
    'frqlim'        'real'       []             [3 srate/2]; ...
    'freqscale'     'string'   {'linear' 'log'}      'log';...
    'pcfac'        'integer'      []             [7];...
    'cycles'      'real'      []             [3 0.5];...
    'timesout'      'integer'      []             [200];...
    'nfreqs'       'integer'      []             [];...
    'padratio'     'integer'      []             [4];...
    'winsize'      'integer'      []             [];...
    'wletmethod'   'string'    {'dftfilt' 'dftfilt2' 'dftfilt3'}   'dftfilt3';...
    'icatype'      'string'    {'amica' 'infomax'}           'infomax'...
    }, 'inputgui');
if isstr(g), error(g); end;

if ~iscell(eegdata)
    eegtemp{1} = eegdata;
else
    eegtemp = eegdata;
end

%% set timewindow to use


TWindow = g.cycles(1)*(1/g.frqlim(1));

if isempty(g.winsize) || g.winsize < TWindow;
    g.winsize = TWindow;
end

g.winsize = g.winsize*srate;



kt = [];
meanpwrCond = [];
complist = g.plotcomps;
CondTrials = [];
CondTimewindows = [];
CondTimevec = [];
Trials = [];



pwr = [];

for ink = 1:length(eegtemp)
    eeg = eegtemp{ink};
    for cp = 1:length(g.plotcomps)
        [ersp, itc, powbase, timevectmp, freqvec, eboot, pboot, tfdata] = newtimef...
            (eeg(g.plotcomps(cp),:,:), frames, tlimits, srate, g.cycles,...
            'freqs', [g.frqlim(1) g.frqlim(end)], ...
            'freqscale', g.freqscale,'timesout', g.timesout,...
            'wletmethod', g.wletmethod, 'baseline',NaN,... %'nfreqs', length([g.frqlim(1):g.frqlim(end)])*nfreqs,...
            'plotitc', 'off', 'plotersp', 'off','nfreqs', g.nfreqs, 'padratio', g.padratio, 'winsize', g.winsize);        
        pwr1 = tfdata(:,:);
        pwr1 = pwr1.*conj(pwr1);
        pwr1 = 10*log10(pwr1);
        pwr_tmp(:,:,cp) = pwr1'; clear pwr1
        fprintf('\nIC %s done...',g.plotcomps(cp));
        meanpwrCond(ink,cp,:) = squeeze(mean(pwr_tmp(:,:,cp),1));
    end;
    ntw_trial = size(tfdata,2);
    n_trials(ink) = size(tfdata,3);
    pwr = cat(1,pwr,pwr_tmp);
    timevec_trial = timevectmp;
    pwr_tmp = [];
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate and remove mean spectrum from all comps   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for cp = 1:size(pwr,3)
    meanpwr(cp,:) = squeeze(mean(pwr(:,:,cp),1));
    pwr(:,:,cp) = pwr(:,:,cp) - repmat(meanpwr(cp,:),[size(pwr,1) 1]); % remove mean
end;

clear alldat

pwr = reshape(pwr,size(pwr,1),size(pwr,2)*size(pwr,3));% windows x comps*freqs

fprintf('\nRemoving mean for each spectral window...\n');
for tt = 1:size(pwr,1)
    rowmeans(tt,:) = mean(pwr(tt,:)); % collect the means taken out for future reference
    pwr(tt,:) = pwr(tt,:) - mean(pwr(tt,:)); % take out mean of each row (windows)
end;


pcs =  round(sqrt(size(pwr,2)/g.pcfac));

%--------------------------------------------
fprintf('\nPCA''ing to %s dimensions.\n',int2str(pcs));%%%%%%%%%%%%%%%%%%%%%%
[U,S,V] = svds(pwr',pcs);% if you scale 'acts', you can't scale the 'eigvec'
pc = (U*S)'; % scale 'activations' for appropriate weighting in decomp of pc
eigvec = V;


numrows = size(pwr,1);
numframes = size(pwr,2);
clear pceigvec


%% run ICA

if strcmp(g.icatype, 'amica')
    
    indPlugin = 0;
    
    try, PLUGINLIST = evalin('base', 'PLUGINLIST'); catch, PLUGINLIST = []; end;
    if ~isempty(PLUGINLIST) && isfield(PLUGINLIST, 'plugin')
        indPlugin = strmatch(lower('amica'), lower({ PLUGINLIST.plugin }), 'exact');
    end
    
    if indPlugin == 0;
        installRes = plugin_askinstall('ICLabel');
        if installRes == 0
            errordlg2('Cannot run AMICA without amica plugin installed. Running infomax ICA instead. To run AMICA please install AMICA plugin using the eeglab extension manager.',...
                'AMICA plugin is not installed');
            
            fprintf('\n Running infomax ICA...\n');
            
            [IMAweights,IMAsphere] = runica(pc,'extended',1,'maxsteps',2000,'stop',1e-8);
            
        end
    end
    
    
    [IMAweights,IMAsphere,mods] = runamica15(pc,...
        'num_chans',size(pc,1),...
        'num_models', 1,...
        'max_iter', 2000,...
        'do_reject',1,...
        'numrej', 5,...
        'rejsig', 5);
    
else
    
    [IMAweights,IMAsphere] = runica(pc,'extended',1,'maxsteps',2000,'stop',1e-8);
    
end



%%%%%

