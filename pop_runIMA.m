%
% selects components using IC label, performs time-frequency decomposition, PCA and runs Independent Modulator
% Analysis (IMA) for a single subject with a single condition 
%
%[EEG] = pop_runIMA(EEGc,varargin)
%
% Author: Johanna Wagner, Swartz Center for Computational Neuroscience, UC San Diego, 2019
% adapted from a function written by Julie Onton
%
%% Example computing IMA a single subjects selecting brain and muscle ICs using IC label, with
% parameters for time-frequency decomposition: log scale, frequency limit 6
% to 120Hz, wavelet cycles [6 0.5], reducing the dimensions of timewindows
% of the tf decomposition using pfac 8 (for a description of pfac see
% below), and using AMICA as the ica algorithm for IMA
% 
%  >>  [EEG] = pop_runIMA_study(EEG, 'selectICs', {'brain' 'muscle'}, 'freqscale', 'log', 'frqlim', [6 120], 'cycles', [6 0.5], 'pcfac', 8,...
%      'icatype', 'amica') 

%
% INPUTS:
% EEGc -- a single EEG structure
% plotcomps -- [vector of IC indices] components to run IMA on (or select via IC label -- see option 'selectICs' below)
% frqlim -- [minfrq maxfrq] minimum and maximum frequencies to include in spectral decomposition
% freqscale -- ['linear',or 'log']
% pcfac -- [integer] - regulates the dimension that the time windows of the spectral data will
%                      be reduced to before IMA - the smaller pcfac, the more dimensions will be
%                      retained ndims = (freqs*ICs)/pcfac where freqs is the number of estimated
%                      frequencies and ICs is the number of ICs (default is 7)
% epochlength -- [number] if data is not epoched size of data epoch in seconds (default 4) 
% selectICs   -- ['cell'] select ICs with IC label can be   {'brain' 'muscle' 'eye' 'artefact' 'off'}
%                         several categories can be chosen at once   default  {'off'};...
% iclthreshold_{brain,musecle, eye, artefact} -- [number between 0 and 1] defines the thersold for IC label to select the components 
%                     specified in 'selectICs' a higher number indicates a
%                     stricter threshold  default [0.7].
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
% EEG structure pointing to the .ima file at EEG.etc.IMA.filename and EEG.etc.IMA.filepath with the IMA structure
% IMA - structure with fields of IMA output (see detailed description below)
% results are saved in the original subject folder with the extension .ima 
% this file has the same properties as a .mat file and can be loaded in Matlab using 
% load(['filename.ima'], '-mat' )
%
%detailed description of IMA outputs:
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
% IMA.winsize - window length for computing spectra in the time-frequency decomposition
% IMA.epochlength - epochlength used for computing time-frequency decomposition in seconds*
% IMA.eigvec - pc backprojection in time
% IMA.pc - pc spectral backprojection
% IMA.timefreq - time-frequency decomposition (spectograms for each IC)
% IMA.timepntCond - total number of timepoints in time-frequency decomposition
% IMA.timevec - timevector of full length of time-frequency decomposition
% IMA.subjfilename - the filename of the .ima file
% IMA.subjfilepath - the filepath of .ima file

function [EEG, IMA] = pop_runIMA(EEGc,varargin)
IMA = []; EEG = EEGc;


g = finputcheck(varargin, { 'plotcomps'     'integer'   []             []; ...
    'frqlim'        'real'       []             [3 EEGc.srate/2]; ...
    'freqscale'     'string'   {'linear' 'log'}      'log';...
    'pcfac'        'integer'      []             [7];...
    'cycles'      'real'      []             [3 0.5];...
    'overlapfac'      'real'      []            [0.6];...
    'nfreqs'       'integer'      []             [];...
    'padratio'     'integer'      []             [4];...
    'wletmethod'   'string'    {'dftfilt' 'dftfilt2' 'dftfilt3'}   'dftfilt3';...
    'winsize'      'integer'      []             [];...
    'epochlength' 'real'          []             [6];...
    'selectICs'   'cell'     {'brain' 'muscle' 'eye' 'artefact' 'off'}     {'off'};...
    'iclthreshold_brain' 'real'    []     [0.7];...
    'iclthreshold_muscle' 'real'    []     [0.7];...
    'iclthreshold_eye' 'real'    []     [0.7];...
    'iclthreshold_artefact' 'real'    []     [0.7];...
    'usedip'       'string'  {'on' 'off'}       'on';...     
    'dipthreshold' 'real'    []     [0.15];...
    'icatype'      'string'    {'amica' 'infomax'}           'infomax'...
    }, 'inputgui');
if isstr(g), error(g); end;

% Check if ICLABEL is installed
try PLUGINLIST = evalin('base', 'PLUGINLIST'); catch, PLUGINLIST = []; end
if ~isempty(PLUGINLIST) && isfield(PLUGINLIST, 'plugin')
    indPlugin = strmatch(lower('ICLabel'), lower({ PLUGINLIST.plugin }), 'exact');
end

flag_iclabel = 1;
if indPlugin == 0
    flag_iclabel = plugin_askinstall('ICLabel');
end


%%  Display GUI
if nargin ==1
    
    freqlim_def = [1 EEGc.srate/2-10];
    iclabel_list = {'Brain', 'Muscle', 'Eye', 'Heart'};  
    freqscale_list = {'linear' 'log'} ;
    
%     'chbx_brain', 'chbx_muscle'. 'chbx_eye', 'chbx_heart', 'ed_brainthrs', 'ed_musclethrs', 'ed_eyethrs', 'ed_heartthrs'

    if flag_iclabel
        cb_chbx_icindx  = ['val_chbx_icindx = get(findobj(''Tag'', ''chbx_icindx''), ''Value'');' ...
            'if val_chbx_icindx,'...
            'set(findobj(''Tag'', ''chbx_iclabel''), ''Value'', 0);'...
            'set(findobj(''Tag'', ''ed_icindx''), ''enable'', ''on'');'...
            'set(findobj(''Tag'', ''chbx_brain''), ''enable'', ''off'');'...
            'set(findobj(''Tag'', ''chbx_eye''), ''enable'', ''off'');'...
            'set(findobj(''Tag'', ''chbx_muscle''), ''enable'', ''off'');'...
            'set(findobj(''Tag'', ''chbx_heart''), ''enable'', ''off'');'...
            'set(findobj(''Tag'', ''ed_brainthrs''), ''enable'', ''off'');'...
            'set(findobj(''Tag'', ''ed_musclethrs''),''enable'', ''off'');'...
            'set(findobj(''Tag'', ''ed_eyethrs''), ''enable'', ''off'');'...
            'set(findobj(''Tag'', ''ed_heartthrs''), ''enable'', ''off'');'...
            'set(findobj(''Tag'', ''chbx_icindx''), ''Value'', 1);'...
            'else;'...
            'set(findobj(''Tag'', ''chbx_iclabel''), ''Value'', 1);'...
            'set(findobj(''Tag'', ''ed_icindx''), ''enable'', ''off'');'...
            'set(findobj(''Tag'', ''chbx_brain''), ''enable'', ''on'');'...
            'set(findobj(''Tag'', ''chbx_muscle''), ''enable'', ''on'');'...
             'set(findobj(''Tag'', ''chbx_eye''), ''enable'', ''on'');'...
            'set(findobj(''Tag'', ''chbx_heart''), ''enable'', ''on'');'...
            'set(findobj(''Tag'', ''ed_brainthrs''), ''enable'', ''on'');'...
            'set(findobj(''Tag'', ''ed_musclethrs''),''enable'', ''on'');'...
            'set(findobj(''Tag'', ''ed_eyethrs''), ''enable'', ''on'');'...
            'set(findobj(''Tag'', ''ed_heartthrs''), ''enable'', ''on'');'...
            'set(findobj(''Tag'', ''chbx_icindx''), ''Value'', 0);'...
            'end;'];
        
        cb_chbx_iclabel = ['val_chbx_iclabel = get(findobj(''Tag'', ''chbx_iclabel''), ''Value'');' ...
            'if ~val_chbx_iclabel,'...
            'set(findobj(''Tag'', ''chbx_iclabel''), ''Value'', 0);'...
            'set(findobj(''Tag'', ''ed_icindx''), ''enable'', ''on'');'...
             'set(findobj(''Tag'', ''chbx_brain''), ''enable'', ''off'');'...
            'set(findobj(''Tag'', ''chbx_muscle''), ''enable'', ''off'');'...
             'set(findobj(''Tag'', ''chbx_eye''), ''enable'', ''off'');'...
            'set(findobj(''Tag'', ''chbx_heart''), ''enable'', ''off'');'...
            'set(findobj(''Tag'', ''ed_brainthrs''), ''enable'', ''off'');'...
            'set(findobj(''Tag'', ''ed_musclethrs''),''enable'', ''off'');'...
            'set(findobj(''Tag'', ''ed_eyethrs''), ''enable'', ''off'');'...
            'set(findobj(''Tag'', ''ed_heartthrs''), ''enable'', ''off'');'...
            'set(findobj(''Tag'', ''chbx_icindx''), ''Value'', 1);'...
            'else;'...
            'set(findobj(''Tag'', ''chbx_iclabel''), ''Value'', 1);'...
            'set(findobj(''Tag'', ''ed_icindx''), ''enable'', ''off'');'...
            'set(findobj(''Tag'', ''chbx_brain''), ''enable'', ''on'');'...
            'set(findobj(''Tag'', ''chbx_muscle''), ''enable'', ''on'');'...
            'set(findobj(''Tag'', ''chbx_eye''), ''enable'', ''on'');'...
            'set(findobj(''Tag'', ''chbx_heart''), ''enable'', ''on'');'...
            'set(findobj(''Tag'', ''ed_brainthrs''), ''enable'', ''on'');'...
            'set(findobj(''Tag'', ''ed_musclethrs''),''enable'', ''on'');'...
            'set(findobj(''Tag'', ''ed_eyethrs''), ''enable'', ''on'');'...
            'set(findobj(''Tag'', ''ed_heartthrs''), ''enable'', ''on'');'...
            'set(findobj(''Tag'', ''chbx_icindx''), ''Value'', 0);'...
            'end;'];
    else
        cb_chbx_icindx = ['set(findobj(''Tag'', ''chbx_iclabel''), ''Value'', 0);'...
            'set(findobj(''Tag'', ''ed_icindx''), ''enable'', ''on'');'...
            'set(findobj(''Tag'', ''chbx_brain''), ''enable'', ''off'');'...
            'set(findobj(''Tag'', ''chbx_muscle''), ''enable'', ''off'');'...
            'set(findobj(''Tag'', ''chbx_eye''), ''enable'', ''off'');'...
            'set(findobj(''Tag'', ''chbx_heart''), ''enable'', ''off'');'...
            'set(findobj(''Tag'', ''ed_brainthrs''), ''enable'', ''off'');'...
            'set(findobj(''Tag'', ''ed_musclethrs''),''enable'', ''off'');'...
            'set(findobj(''Tag'', ''ed_eyethrs''), ''enable'', ''off'');'...
            'set(findobj(''Tag'', ''ed_heartthrs''), ''enable'', ''off'');'...
            'set(findobj(''Tag'', ''chbx_icindx''), ''Value'', 1);'];
        
        cb_chbx_iclabel = ['set(findobj(''Tag'', ''chbx_iclabel''), ''Value'', 0);'...
            'set(findobj(''Tag'', ''ed_icindx''), ''enable'', ''on'');'...
            'set(findobj(''Tag'', ''chbx_brain''), ''enable'', ''on'');'...
            'set(findobj(''Tag'', ''chbx_muscle''), ''enable'', ''on'');'...
            'set(findobj(''Tag'', ''chbx_eye''), ''enable'', ''on'');'...
            'set(findobj(''Tag'', ''chbx_heart''), ''enable'', ''on'');'...
            'set(findobj(''Tag'', ''ed_brainthrs''), ''enable'', ''on'');'...
            'set(findobj(''Tag'', ''ed_musclethrs''),''enable'', ''on'');'...
            'set(findobj(''Tag'', ''ed_eyethrs''), ''enable'', ''on'');'...
            'set(findobj(''Tag'', ''ed_heartthrs''), ''enable'', ''on'');'...
            'set(findobj(''Tag'', ''chbx_icindx''), ''Value'', 1);'];
    end
    
    
    
    %%
    uilist = {{'style' 'text' 'string' 'Select ICs'} ...
        { 'style' 'Checkbox'   'string' 'IC indices' 'Tag' 'chbx_icindx', 'callback', cb_chbx_icindx 'Value' 1} {'style'  'edit' 'string' ' ' 'tag' 'ed_icindx'}...
        { 'style' 'Checkbox'   'string' 'ICLabel tags' 'Tag', 'chbx_iclabel', 'callback', cb_chbx_iclabel 'Value' 0} {'style' 'text' 'string' 'Labels'} {'style' 'text' 'string' 'Threshold (%)'}...
        { 'style' 'Checkbox'   'string' iclabel_list{1} 'Tag', 'chbx_brain', 'Value' 0,  'enable', 'off'}   {'style'  'edit' 'string' num2str(g.iclthreshold_brain) 'tag' 'ed_brainthrs',  'enable', 'off'}...
        { 'style' 'Checkbox'   'string' iclabel_list{2} 'Tag', 'chbx_muscle', 'Value' 0,  'enable', 'off'}  {'style'  'edit' 'string' num2str(g.iclthreshold_muscle) 'tag' 'ed_musclethrs' 'enable', 'off'}...
        { 'style' 'Checkbox'   'string' iclabel_list{3} 'Tag', 'chbx_eye', 'Value' 0,  'enable', 'off'}     {'style'  'edit' 'string' num2str(g.iclthreshold_eye) 'tag' 'ed_eyethrs' 'enable', 'off'}...
        { 'style' 'Checkbox'   'string' iclabel_list{4} 'Tag', 'chbx_heart', 'Value' 0,  'enable', 'off'}   {'style'  'edit' 'string' num2str(g.iclthreshold_artefact) 'tag' 'ed_heartthrs' 'enable', 'off'}...
        {}...
        {'style' 'text' 'string' 'Freq. limits (Hz)'} ...
        {'style' 'edit' 'string' num2str(freqlim_def) 'tag' 'freqlimits'}...
        {'style' 'text' 'string' 'Freq. scale'}...
        {'style' 'popupmenu'  'string' freqscale_list 'tag' 'freqscale' 'value' 2}...
        {'style' 'text' 'string' 'pcfac (see Help)'}...
        {'style' 'edit'  'string' num2str(g.pcfac) 'tag' 'pcfac'}...
        {'style' 'text' 'string' 'pop_runima options (see Help)'}...
        {'style' 'edit' 'string' ' ' 'tag' 'ed_opt'} {}};
    
    ht = 10; wt = 3 ;
   geom = {{wt ht [0 0]  [1 1]} ...
          {wt ht [0.1 1]  [1 1]} {wt ht [1 1] [2 1]} ...
          {wt ht [0.1 2]  [1 1]} {wt ht [1 2] [1 1]} {wt ht [1.45 2] [1 1]}...
          {wt ht [1 3] [0.5 1]} {wt ht [1.5 3] [0.3 1]}...
          {wt ht [1 4] [0.5 1]} {wt ht [1.5 4] [0.3 1]}...
          {wt ht [1 5] [0.5 1]} {wt ht [1.5 5] [0.3 1]}...
          {wt ht [1 6] [0.5 1]} {wt ht [1.5 6] [0.3 1]}...
          {wt ht [0.1 7]  [1 1]}...
          {wt ht [0 8]  [1 1]}  {wt ht [0.5 8] [0.5 1]}  {wt ht [1.1 8]  [1 1]}  {wt ht [1.5 8]  [0.5 1]} {wt ht [2 8]  [1 1]}  {wt ht [2.5 8]  [0.5 1]} ...
          {wt ht [0 9]  [1 1]}  {wt ht [1 9] [2 1]} ...
          {wt ht [0.1 10]  [1 1]}...
          };
    
    [result, ~, ~, resstruct, ~] = inputgui('title','Run Independent Modulator Analysis -- pop_runima', 'geom', geom, 'uilist',uilist, 'helpcom','pophelp(''pop_runIMA'');');
    %%
    if isempty(result), return; end;
    if resstruct.chbx_icindx
        if isempty(str2num(resstruct.ed_icindx))
            g.plotcomps = 1:size(EEG.icawinv,2);
        else
            g.plotcomps = str2num(resstruct.ed_icindx);
        end
    else
        g.plotcomps = [];
        
        g.iclthreshold_brain = str2num(resstruct.ed_brainthrs);
        g.iclthreshold_eye = str2num(resstruct.ed_eyethrs);
        g.iclthreshold_muscle = str2num(resstruct.ed_musclethrs);
        g.iclthreshold_artefact = str2num(resstruct.ed_heartthrs);  
    end
    
    if resstruct.chbx_iclabel
        guilabels = [resstruct.chbx_brain resstruct.chbx_muscle resstruct.chbx_eye resstruct.chbx_heart];
        g.selectICs = iclabel_list{logical(guilabels)};
    else
        g.selectICs = 'off';
    end
    
    g.frqlim = str2num(resstruct.freqlimits);
    g.freqscale = freqscale_list{resstruct.freqscale};
    g.pcfac = str2num(resstruct.pcfac);
    
 %%   
    % Retrieve optional parameters
    tmpoptparams   = eval( [ '{' get(findobj(gcf,'tag','ed_opt'),'string') '}' ] );
    tmpparams_name = tmpoptparams(1:2:end);
    
    % Update parameters here
    c =1;
    for i = 1: length(tmpparams_name)
        g.(tmpparams_name{i}) =  tmpoptparams{c+1};
        c = c+2;
    end
    
end

%% selecting ICs with IClabel

    
%% check if IC label plugin is installed
%% ask if plugin should be downloaded
    
if sum(ismember(g.selectICs,'off')) == 0
%     indPlugin = 0;
%     try, PLUGINLIST = evalin('base', 'PLUGINLIST'); catch, PLUGINLIST = []; end;
%     if ~isempty(PLUGINLIST) && isfield(PLUGINLIST, 'plugin')
%         indPlugin = strmatch(lower('ICLabel'), lower({ PLUGINLIST.plugin }), 'exact');
%     end
%     
%     if indPlugin == 0;
%         installRes = plugin_askinstall('ICLabel');
%         if installRes == 0
%             errordlg2('Cannot select ICs with IC Label. Please install ICLabel in previous step or using the eeglab extension manager.',...
%                 'ICLabel not installed');           
%             error('ICLabel not installed. Cannot select ICs with IC Label. Please install ICLabel in previous step or using the eeglab extension manager.')
%         end
%     end
    
    
    fprintf('selecting ICs with IClabel... \n')
    
    %% check if dipoles are computed
    
    if sum(strcmp(g.selectICs,'brain')) == 1 && strcmp(g.usedip, 'on') && ~isfield(EEGc.dipfit,'model');
          errordlg2('You chose to use dipole residual variance in addition to IC Label to select brain components. You need to first compute dipole locations.',...
                'No dipole locations present');           
            error('No dipole locations present. You need to first compute dipole locations.')
    end
    

    %% select ICs with IC label
    
    EEGc = pop_iclabel(EEGc, 'Default'); % run IC label
    
    [VAL,IND] = max(EEGc.etc.ic_classification.ICLabel.classifications(:,1:7),[],2); % check which categories IC have highest classification accuracy
    
    % select braincomps
    if strcmp(g.usedip, 'on') % use IC label classification plus dipole RV for selecting ICs
    LOGIC{1} = find((IND == 1) & (VAL > g.iclthreshold_brain) & [EEGc.dipfit.model.rv]' < g.dipthreshold);
    else
    LOGIC{1} = find((IND == 1) & (VAL > g.iclthreshold_brain)); % use only IC label classification
    end
    
    LOGIC{2} = find((IND == 2) & (VAL > g.iclthreshold_muscle)); % find muscle ICs with classification accuracy above threshold
    LOGIC{3} = find((IND == 3) & (VAL > g.iclthreshold_eye)); % find eye ICs with classification accuracy above threshold
    LOGIC{4} = find((IND == 4) & (VAL > g.iclthreshold_artefact)); % find heart ICs with classification accuracy above threshold
    indICs = find(ismember(lower(g.selectICs), {'brain' 'muscle' 'eye' 'heart'}));
    
    INDKEEP = [];
    for iki = 1:length(indICs);
        INDKEEP = [INDKEEP LOGIC{iki}']; % find which ICs to keep according to predefined categories 
    end
    g.plotcomps = INDKEEP;

end


if isempty(g.plotcomps);
    g.plotcomps = 1:size(EEGc.icaact,1); % check if plotcomps is empty, otherwise use all
end

%% if data is not epoched - estimate overlap of epochs

if length(size(EEGc.data)) == 2;
    
    fprintf('data is not epoched \n')
    fprintf('estimating epoch overlap... \n')
    
     % calculate windowsize needed to estimate n cycles of lowest frequency
     TWindow = (g.cycles(1)*(1/g.frqlim(1)));
    
   if isempty(g.winsize) || g.winsize < TWindow; % check if g.winsize is empty or if predefined winsize is smaller than required winsize for lowest frequency
        g.winsize = TWindow;
   end
   
      overlap = (g.winsize*g.overlapfac); % estimate overlap (in seconds)
      timesout = round(g.epochlength/(g.winsize-overlap)); % number of timewindows fitting in an epoch with given overlap
      overlapframes = overlap*EEGc.srate; % estimate overlap in frames
      
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
         EEGcurr = EEGc;
        
        % to save original timevector of data and track data that will be removed during epoching
        % due to boundaries add data channel with timevector
        
        EEGcurr.data(end+1,:) = EEGcurr.times;
        EEGcurr.nbchan = size(EEGcurr.data,1);
        
        if ~isempty(EEGcurr.chanlocs)
            EEGcurr.chanlocs(end+1).label = 'timevec';
        end;
        
        timevecorig = EEGcurr.times; % save original timevector in variable
        
        % estimate how many epochs fit into the data given predefined epoch length and estimated overlap
        NEO = ceil(EEGcurr.pnts/((g.epochlength*EEGcurr.srate-overlapframes)));
        % estimate timing of markers to add to data for epoching
        LATEO = (overlapframes:g.epochlength*(EEGcurr.srate)-overlapframes:(g.epochlength*(EEGcurr.srate)-overlapframes)*NEO)
        markersEO = [ones(size(LATEO))'*100 LATEO']; % construct marker structure to feed into eeglab
        
        [EEGcurr] = pop_importevent(EEGcurr, 'event', ... % import event markers to eeglab
            markersEO, 'fields', {'type', 'latency'}, ...
            'append', 'yes', 'align', 0, 'timeunit', NaN );
        EEGcurr = eeg_checkset( EEGcurr );
        
        fprintf('epoching data... \n')
        
        EEGcurr = pop_epoch( EEGcurr, {  num2str(100) }, [0  g.epochlength], 'epochinfo', 'yes'); % epoch EEG data
        EEGcurr = eeg_checkset( EEGcurr );
        eegdata{1} = EEGcurr.icaact;             % save IC activities of dataset in cell array
        timevecproc =  EEGcurr.data(end,:,:);    % save processed version of original timevector in variable
        EEGcurr = pop_select( EEGcurr,'nochannel',size(EEGcurr.data,1)); % remove previously added data channel with timevector from EEG structure
         
end


%% run IMA analysis    

[IMAweights,IMAsphere, meanpwr, freqvec, n_trials, ntw_trials,...
    pcs, eigvec, pc, timefreq]  = runIMA(eegdata,EEGcurr.pnts, [EEGcurr.xmin EEGcurr.xmax]*1000, EEGcurr.srate,...
    'plotcomps', g.plotcomps, 'cycles', g.cycles, 'timesout', timesout, 'frqlim', g.frqlim, 'freqscale', g.freqscale, 'winsize', g.winsize, 'icatype', g.icatype);


%%%%%

% resample processed timevector (after epoching) to timepoint output of
% timefrequency runIMA results
timevec_res = downsample(squeeze(timevecproc),size(squeeze(timevecproc),1)/ntw_trials);
timevec_res = timevec_res(:);

%% warp timef data  and eigvec to start and beginning of processed timevector - interpolates missing data

timefreqtmp = [];
eigvectmp = [];
[timefreqtmp] = resample(timefreq,timevec_res);
[eigvectmp] = resample(eigvec,timevec_res);
timepntcond = [1:n_trials*ntw_trials]; % select trials for condition
timevec = timevec_res;


%% save results in IMA structure

IMA.wts = IMAweights;
IMA.sph = IMAsphere;
IMA.meanpwr = meanpwr;
IMA.timevec = timevec;
IMA.freqvec = freqvec;
IMA.freqscale = g.freqscale;
IMA.freqlim = g.frqlim;
IMA.npcs = pcs;
IMA.complist = g.plotcomps;
IMA.srate = EEGc.srate;
IMA.ntrials = n_trials;
IMA.ntw_trials = ntw_trials;
IMA.winsize = g.winsize;
IMA.epochlength = g.epochlength;
IMA.eigvec = eigvectmp;
IMA.pc = pc;
IMA.timefreq = timefreqtmp;
IMA.meanpwrCond = [];
IMA.timepntCond = timepntcond;
IMA.timevec = timevec;
IMA.condition = [];
IMA.subjfilename = {EEGc.filename};
IMA.subjfilepath = {EEGc.filepath};


% save IMA data structure with same name and in same folder as EEG dataset
filename = extractBefore(EEGc.filename,'.set')

save([EEGc.filepath '/' filename '.ima'],'IMA')


% save filename and filepath of IMA data in EEG structure
EEGc.etc.IMA.filename = [filename '.ima'];
EEGc.etc.IMA.filepath = EEGc.filepath;

EEGc = pop_saveset( EEGc, 'filename', [EEGc.filename], 'filepath', EEGc.filepath);    
EEG = EEGc;


fprintf('\ndone.\n');


