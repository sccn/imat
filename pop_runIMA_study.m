%
% Selects components using IC label, performs time-frequency decomposition, PCA and runs Independent Modulator
% Analysis (IMA) over conditions for each subject in a study
%
%[STUDY] = pop_runIMA_study(STUDY,varargin)
%
% Author: Johanna Wagner, Swartz Center for Computational Neuroscience, UC San Diego, 2019
% adapted from a function written by Julie Onton
%
%% Example computing IMA over all subjects contained in study selecting only brain ICs using IC label, with
% parameters for time-frequency decomposition: log scale, frequency limit 6
% to 120Hz, wavelet cycles [6 0.5], reducing the dimensions of timewindows
% of the tf decomposition using pfac 8 (for a description of pfac see
% below), and using AMICA as the ica algorithm for IMA
%
%  >>  [IMA] = pop_runIMA_study(STUDY, 'selectICs', {'brain'}, 'freqscale', 'log', 'frqlim', [6 120], 'cycles', [6 0.5], 'pcfac', 8,...
%      'icatype', 'amica')
%
%
% INPUTS:
% STUDY - STUDY structure with information on stored IMA files
% subject - subject to compute IMA on - provide subject code/number as string i.e. '3'
%           if empty computes results for all subjects in the study
% frqlim -- [minfrq maxfrq] minimum and maximum frequencies to include in spectral decomposition
% freqscale -- ['linear',or 'log']
% pcfac -- [integer] - regulates the dimension that the time windows of the spectral data will
%                      be reduced to before IMA - the smaller pcfac, the more dimensions will be
%                      retained ndims = (freqs*ICs)/pcfac where freqs is the number of estimated
%                      frequencies and ICs is the number of ICs (default is 10)
% epochlength -- [number] if data is not epoched size of data epoch in seconds (default 4)
% selectICs   -- ['cell'] select ICs with IC label can be   {'brain' 'muscle' 'eye' 'artefact' 'off'}
%                         several categories can be chosen at once   default  {'off'};...
% iclthreshold -- [number between 0 and 1] defines the thersold for IC label to select the components
%                     specified in 'selectICs' a higher number indicates a stricter threshold  default [0.7]
% usedip      -- ['string'] {'on' 'off'} use dipole residual variance (RV) in addition to ICLabel to select
%                           brain components, default {'on'}
% dipthreshold -- 'real'    [number between 0 and 1] which cutoff to use for dipole RV for selecting brain
%                           components default [0.15];...
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
% overlapfac -- [number between 0 and 1] percentage by which to overlap
%                       timewindows a larger number means larger overlap - default 0.6
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
%                       for more details {default: 'dftfilt3'}
% winsize  --   [integer] in sec size of timwindows to use for time-frequency decomposition
%                         needs to be smaller than epochsize. Note: this parameter
%                       is overwritten when the minimum frequency and number of cycles requires
%                       a longer time window. {default: cycles*(length of a cycle at lowest freq in sec)}
%
%
%%% OUTPUT
% STUDY structure with fields under STUDY.etc.IMA.imafilename and STUDY.etc.IMA.imafilepath
% pointing to the .ima files and filepaths for each subject;
% IMA - structure with fields of IMA output (see detailed description below)
% results are saved subjectwise in the original subject folder with the extension .ima
% this file has the same properties as a .mat file and can be loaded in Matlab using
% load(['filename.ima'], '-mat' )
%
%
% Detailed description of IMA output structure:
% IMA.wts - weights of IMA decomposition
% IMA.sph - spheres of IMA decomposition
% IMA.meanpwr - mean power spectra of single ICs
% IMA.timevec - timevector
% IMA.freqvec - frequency vector
% IMA.freqscale - frequency scale of computed spectra ('log' or 'linear')
% IMA.freqlim - frequency limits of spectra
% IMA.npcs - number of dimensions that have been used to reduce the data before IMA
% IMA.complist - component indices on which IMA was run on
% IMA.srate - original sampling rate of EEG data used to compute the spectra
% IMA.ntrials - number of trials used to commpute the tie-frequency decomposition
% IMA.ntw_trials - number of timewindows per trial
% IMA.eigvec - pc backprojection in time
% IMA.pc - pc spectral backprojection
% IMA.timefreq - time-frequency decomposition (spectograms for each IC)
% IMA.timepntCond - total number of timepoints in time-frequency
%                   decomposition in each condition
% IMA.timevec_cond - timevector of full length of time-frequency
%                     decomposition for each condition
% IMA.meanpwrCond - mean power spectra for each IC and each condition
% IMA.condition - names and order of conditions
% IMA.STUDYname - filename of the STUDY the IMA decomposition belongs to
% IMA.STUDYfilepath - filepath of the STUDY the IMA decomposition belongs to
% IMA.subj - subject the IMA has been comuted on
% IMA.subjfilename - filenames of the EEG data the IMA has been computed on
% IMA.subjfilepath - filepath of the EEG data the IMA has been computed on


function [STUDY] = pop_runIMA_study(STUDY,ALLEEG, varargin)


g = finputcheck(varargin, {'frqlim'        'real'       []             []; ...
    'freqscale'     'string'   {'linear' 'log'}      'log';...
    'pcfac'        'integer'      []             [];...
    'cycles'      'real'      []             [3 0.5];...
    'overlapfac'      'real'      []            [0.6];...
    'nfreqs'       'integer'      []             [];...
    'padratio'     'integer'      []             [4];...
    'wletmethod'   'string'    {'dftfilt' 'dftfilt2' 'dftfilt3'}   'dftfilt3';...
    'winsize'      'integer'      []             [];...
    'epochlength' 'real'          []             [6];...
    'selectICs'   'cell'     {'brain' 'muscle' 'eye' 'artefact' 'off'}     {'off'};...
    'iclthreshold' 'real'    []     [0.7];...
    'usedip'       'string'  {'on' 'off'}       'on';...
    'dipthreshold' 'real'    []     [0.15];...
    'icatype'      'string'    {'amica' 'infomax'}           'infomax'...
    }, 'inputgui');
if isstr(g), error(g); end;

% if narging==1
%     %Display GUI
%
% end


%% check format of EEG set
ncond = length({STUDY.datasetinfo.subject})/length(unique({STUDY.datasetinfo.subject}));
nsubj = length(STUDY.subject);

for iko = 1:nsubj
    indsj = find(ismember({STUDY.datasetinfo.subject}, STUDY(iko).subject));
    for ika = indsj;
        EEGtmp{ika} = pop_loadset('filename',STUDY.datasetinfo(ika).filename,'filepath',STUDY.datasetinfo(ika).filepath);
    end
    
    %% check if IC label plugin is installed
    %% ask if plugin should be downloaded
    
    
    g.plotcomps = [];
    if sum(ismember(g.selectICs,'off')) == 0;
        indPlugin = 0;
        
        try, PLUGINLIST = evalin('base', 'PLUGINLIST'); catch, PLUGINLIST = []; end;
        if ~isempty(PLUGINLIST) && isfield(PLUGINLIST, 'plugin')
            indPlugin = strmatch(lower('ICLabel'), lower({ PLUGINLIST.plugin }), 'exact');
        end
        
        if indPlugin == 0;
            installRes = plugin_askinstall('ICLabel');
            if installRes == 0
                errordlg2('Cannot select ICs with IC Label. Please install ICLabel in previous step or using the eeglab extension manager.',...
                    'ICLabel not installed');
                error('ICLabel not installed. Cannot select ICs with IC Label. Please install ICLabel in previous step or using the eeglab extension manager.')
            end
        end
        
        %% check if dipoles are computed
        fprintf('selecting ICs with IClabel... \n')
        
        if sum(strcmp(g.selectICs,'brain')) == 1 && strcmp(g.usedip, 'on') && ~isfield(EEGtmp{1}.dipfit,'model');
            errordlg2('You chose to use dipole residual variance in addition to IC Label to select brain components. You need to first compute dipole locations.',...
                'No dipole locations present');
            error('No dipole locations present. You need to first compute dipole locations.')
        end
        
        
        %% run IC label
        
        EEGtmp{1} = pop_iclabel( EEGtmp{1}, 'Default');
        
        % check which IC label is most probable for each IC
        [VAL,IND] = max(EEGtmp{1}.etc.ic_classification.ICLabel.classifications(:,1:7),[],2);
        
        % identify braincomponents
        if strcmp(g.usedip, 'on') % based on dipole RV and IC label classification accuracy
            LOGIC{1} = find((IND == 1) & (VAL > g.iclthreshold) & [EEGtmp{1}.dipfit.model.rv]' < g.dipthreshold);
        else
            LOGIC{1} = find((IND == 1) & (VAL > g.iclthreshold)); % based on IC label classification accuracy
        end
        
        
        LOGIC{2} = find((IND == 2) & (VAL > g.iclthreshold)); % identify muscle ICS based on IC Label
        LOGIC{3} = find((IND == 3) & (VAL > g.iclthreshold)); % identify eye ICS based on IC Label
        LOGIC{4} = find((IND == 4) & (VAL > g.iclthreshold)); % identify heart ICS based on IC Label
        indICs = find(ismember(g.selectICs, {'brain' 'muscle' 'eye' 'heart'}));
        
        INDKEEP = [];
        for iki = 1:length(indICs)
            INDKEEP = [INDKEEP LOGIC{iki}']; % check which ICs are part of the predefined categories
        end
        g.plotcomps = INDKEEP;
        
        % update component indices in the STUDY structure
        for indd = indsj;
            STUDY.datasetinfo(indd).comps = INDKEEP;
        end
    end
    
    if isempty(g.plotcomps);
        g.plotcomps = STUDY.datasetinfo(ika).comps;
    end
    
    if isempty(g.frqlim)
        g.frqlim = [3 EEGtmp{1}.srate/2];
    end
    
    %% if data is not epoched - estimate overlap of epochs
    timevecorig = [];
    if length(size(EEGtmp{1}.data)) == 2;
        
        fprintf('data is not epoched \n')
        fprintf('estimating epoch overlap... \n')
        
        % calculate windowsize needed to estimate n cycles of lowest frequency
        TWindow = (g.cycles(1)*(1/g.frqlim(1)));
        
        if isempty(g.winsize) || g.winsize < TWindow;
            g.winsize = TWindow;
        end
        
        overlap = (g.winsize*g.overlapfac);
        timesout = round(g.epochlength/(g.winsize-overlap));
        overlapframes = overlap*EEGtmp{1}.srate;
        
        %% check if windowsize exceeds epoch length or if epoch length is large enough to estimate lowest frequency with n cycles
        if g.winsize >= g.epochlength && g.winsize > TWindow;
            errordlg2('Window size selected for time-frequency decomposition is larger than epoch size. Either choose a smaller time window size or a larger epoch size.',...
                'Window size exceeds epoch limits');
            error('Window size exceeds epoch limits. Window size selected for time-frequency decomposition is larger than epoch size. Either choose a smaller time window size or a larger epoch size.')
        elseif g.winsize >= g.epochlength && g.winsize == TWindow;
            errordlg2('Window size determined for for estimating lowest frequency with n cycles is larger than epoch size. You need to choose a larger epoch size.',...
                'Epoch size is too small for estimating lowest frequency with n cycles');
            error('poch size is too small for estimating lowest frequency with n cycles. Window size determined for for estimating lowest frequency with n cycles is larger than epoch size. You need to choose a larger epoch size.')
        end
        
        
        %% epoching data
        
        for inc = 1:length(EEGtmp)
            EEGcurr = EEGtmp{inc};
            
            % add data channel with timevector to track later rejected parts
            % during epoching
            
            EEGcurr.data(end+1,:) = EEGcurr.times;
            EEGcurr.nbchan = size(EEGcurr.data,1);
            if ~isempty(EEGcurr.chanlocs)
                EEGcurr.chanlocs(end+1).label = 'timevec';
            end;
            
            timevecorig{inc} = EEGcurr.times;
            
            % estimate how may epochs of predefined length fit into data given
            % predefined overlap
            NEO = ceil(EEGcurr.pnts/((g.epochlength*EEGcurr.srate-overlapframes)));
            % determine latency of events for epochs
            LATEO = (overlapframes:g.epochlength*(EEGcurr.srate)-overlapframes:(g.epochlength*(EEGcurr.srate)-overlapframes)*NEO)
            % create event structure with type and latency
            markersEO = [ones(size(LATEO))'*inc*100 LATEO'];
            
            [EEGcurr] = pop_importevent(EEGcurr, 'event', ... % import events to eeg dataset
                markersEO, 'fields', {'type', 'latency'}, ...
                'append', 'yes', 'align', 0, 'timeunit', NaN );
            EEGcurr = eeg_checkset( EEGcurr );
            
            fprintf('epoching data... \n')
            % epoch dataset
            EEGcurr = pop_epoch( EEGcurr, {  num2str(inc*100) }, [0  g.epochlength], 'epochinfo', 'yes');
            EEGcurr = eeg_checkset( EEGcurr );
            eegdata{inc} = EEGcurr.icaact; % save icaactivity of dataset in cell array
            timevecproc{inc} =  EEGcurr.data(end,:,:); % save timevector of processed epoched data to cell array
            EEGcurr = pop_select( EEGcurr,'nochannel',size(EEGcurr.data,1)); %delete perviously created timevector channel
            EEGtemp{inc} = EEGcurr;
            
        end
    else   %% in case that EEG datasets are already epoched EEG structure and icaactivity of dataset in cell array
        for inc = 1:length(EEGtmp)
            EEGtemp{inc} = EEGtmp{inc};
            eegdata{inc} = EEGtmp{inc}.icaact;
        end
    end
    
    %% run IMA
    
    [IMAweights,IMAsphere, meanpwr, freqvec, n_trials, ntw_trials,...
        pcs, eigvec, pc, timefreq,...
        meanpwrCond] = runIMA(eegdata,EEGtemp{1}.pnts, [EEGtemp{1}.xmin EEGtemp{1}.xmax]*1000, EEGtemp{1}.srate,...
        'plotcomps', g.plotcomps, 'cycles', g.cycles, 'timesout', timesout, 'frqlim', g.frqlim, 'freqscale', g.freqscale, 'winsize', g.winsize, 'icatype', g.icatype);
    
    
    %%%%%  downsample the previously saved timevector of the processed dataset to the number of time frequency output timepoints of IMA
    timevec_res = [];
    if length(size(EEGtmp{1}.data)) == 2;
        for inc = 1:length(EEGtmp)
            timevec_res{inc} = downsample(squeeze(timevecproc{inc}(1,1:end-overlapframes,:)),round(size(squeeze(timevecproc{inc}(1,1:end-overlapframes,:)),1)/ntw_trials));
        end
        
        % for inc = 1:length(EEGtmp)
        % timevec_res{inc} = downsample(squeeze(timevecorig{inc}),n_trials(inc)*ntw_trials);
        % end
        
        condvec = [1 n_trials*ntw_trials];
        
        %% warp timef data to original datalength
        
        timefreqtmp = [];
        eigvectmp = [];
        timevec = [];
        opi = 0;
        for inc = 1:length(EEGtmp)
            [timefreqtmp] = [timefreqtmp; resample(timefreq((condvec(inc)+1):(condvec(inc+1)+condvec(inc)),:),timevec_res{inc}(:))];
            [eigvectmp] = [eigvectmp; resample(eigvec((condvec(inc)+1):(condvec(inc+1)+condvec(inc)),:),timevec_res{inc}(:))];
            timepntcond{inc} = (condvec(inc)+1):(condvec(inc+1)+condvec(inc)); % select trials for condition
            timevec = [timevec; opi + timevec_res{inc}(:)];
            opi = timevec_res{inc}(end);
        end
        timefreq = timefreqtmp;
        eigvec = eigvectmp;
    end
    
    
    %% save IMA results in IMA structure
    
    IMA.wts = IMAweights;
    IMA.sph = IMAsphere;
    IMA.meanpwr = meanpwr;
    IMA.freqvec = freqvec;
    IMA.timevec = timevec;
    IMA.timevec_cond = timevec_res;
    IMA.freqscale = g.freqscale;
    IMA.freqlim = g.frqlim;
    IMA.npcs = pcs;
    IMA.complist = g.plotcomps;
    IMA.srate = EEGtemp{1}.srate;
    IMA.ntrials = n_trials;
    IMA.ntw_trials = ntw_trials;
    IMA.eigvec = eigvectmp;
    IMA.pc = pc;
    IMA.timefreq = timefreqtmp;
    IMA.meanpwrCond = meanpwrCond;
    IMA.timepntCond = timepntcond;
    IMA.condition = STUDY.condition;
    IMA.STUDYname = STUDY.filename;
    IMA.STUDYfilepath = STUDY.filepath;
    IMA.subj = STUDY(iko).subject;
    IMA.subjfilename = {STUDY.datasetinfo(indsj).filename};
    IMA.subjfilepath = {STUDY.datasetinfo(indsj).filepath};

    newfilename = extractBefore(STUDY.filename,'.study');
   
    IMA.filename = [IMA.subj{1} '_' newfilename '.ima'];
    
    %% save IMA results to file
    save([EEGtmp{1}.filepath '/' IMA.subj{1} '_' newfilename '.ima'],'IMA');
    
    STUDY.etc.IMA.subjfilename(iko,:) = IMA.subjfilename;
    STUDY.etc.IMA.subjfilepath(iko,:) = IMA.subjfilepath;
    STUDY.etc.IMA.imafilename(iko,:) = [IMA.subj{1} '_' newfilename '.ima'];
    STUDY.etc.IMA.imafilepath(iko,:) = IMA.subjfilepath{1};
    
    
    [STUDY] = pop_savestudy( STUDY, ALLEEG, 'savemode','resave');
     
    clear EEGtmp EEGcurr
    
    
end
fprintf('\ndone.\n');


