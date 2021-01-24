
% collects spectral templates for clustering removing templates whose spectral amplitude are smaller than 
% 1/3 of the maximum template amplitude. Templates are normalized by rms.
% Allows to select templates that are active/ have peaks in a certain range
% of frequencies and warps these spectra to a common target peak over
% subjects
%
% [IMA] = pop_collecttemplates(STUDY,varargin);
%
%
% Author: Johanna Wagner, Swartz Center for Computational Neuroscience, UC San Diego, 2019
% adapted from a function written by Julie Onton
%
%
% Example: chooses templates which have a peak in the range [8 14] and
% stretches them to the target frequency 10 using the subject specific median spectra peak in the range [8 14] 
% >> [IMA] = pop_collecttemplates(STUDY, 'peakfreq', [8 14],'stretch_spectra', 'on', 'targetpeakfreq', 10, 'plot_templ', 'on');
%
% Example: chooses all templates whose spectral amplitude are larger than 
% 1/3 of the maximum template amplitude
% >> [IMA] = pop_collecttemplates(STUDY);
%
%
%
% INPUTS:
% required Inputs:
% STUDY - STUDY structure with information on stored IMA files
%
% optional Inputs:
% peakrange -  'vector of two integers' [min max] include only templates that are active/ have peaks 
%               in this range of frequencies. If empty choose templates with activations in any 
%               frequency range 
% stretch_spectra -  [string] {'on' 'off'} if 'on' stretches templates in
%                   the frequency range defined by 'peakrange' using the
%                   subject specific median template peak frequency in the 'peakrange' to stretch
%                   spectra to a predefined peak as defined in 'targetpeakfreq'.
%                   default is 'off'. Use is only reccomended when a narrow
%                   enough 'peakrange' is defined. i.e. 8-14Hz                    
% targetpeakfreq  -  ['integer'] defines the target frequency to stretch spectra to
%                     when 'stretch_spectra' is 'on', if 'targetpeakfreq' is empty uses
%                    the center frequency of 'peakrange'
% plot_templ - [string] {'on' 'off'} plots overall subject IM templates, the frequency range 
%               selected (if 'peakrange') is given and the selected templates 


function pop_collecttemplates(STUDY,varargin);


g = finputcheck(varargin, {'peakrange'        'integer'       []             []; ...
    'stretch_spectra'     'string'     {'on' 'off'}         'off';...
    'targetpeakfreq'      'integer'       []             []; ...
    'plot_templ'             'string'    {'on' 'off'}    'off';...
    }, 'inputgui');
if isstr(g), error(g); end;

if nargin ==1
     %streatchval = {'on', 'off'}; Check this for removal
     cb_strecht = 'cboxval = get(findobj(''Tag'', ''chbx_strectspec''), ''value''); if cboxval, set(findobj(''Tag'', ''ed_peakfreq''), ''enable'', ''on''); else, set(findobj(''Tag'', ''ed_peakfreq''), ''enable'', ''off'');end;';
     uilist = {{'style' 'text' 'string' 'Freq. range (Hz)'} {'style' 'edit'  'string' ' ' 'tag' 'ed_freqrange'}...
              {'style' 'text' 'string' 'Warp spectra'} {'style' 'checkbox'  'value' 0 'tag' 'chbx_strectspec' 'Callback' cb_strecht}...
              {'style' 'text' 'string' 'Target peak freq. (Hz)'} {'style' 'edit'  'string' ' ' 'Tag' 'ed_peakfreq' 'enable', 'off'}};
    
    ht = 3; wt = 2 ;
    geom = {{wt ht [0 0]  [1 1]} {wt ht [1 0]  [1 1]}...
            {wt ht [0 1]  [1 1]} {wt ht [1 1]  [1 1]}...
            {wt ht [0 2]  [1 1]} {wt ht [1 2]  [1 1]}};
   
    [result, ~, ~, resstruct, ~] = inputgui('title','Collect templates for clustering -- pop_collecttemplates', 'geom', geom, 'uilist',uilist, 'helpcom','pophelp(''pop_collecttemplates'');');
    if isempty(result), return; end;
    g.peakrange = str2num(resstruct.ed_freqrange);
    if resstruct.chbx_strectspec
        g.stretch_spectra = 'on';
    else
        g.stretch_spectra  = 'off'
    end
    g.targetpeakfreq = str2num(resstruct.ed_peakfreq);
end


subjcode = STUDY.subject; 


for iko = 1:length(subjcode)
indsj = find(ismember({STUDY.datasetinfo.subject}, STUDY(iko).subject));

%% load IMA file for current subject
load([STUDY.etc.IMA.imafilepath{iko} filesep STUDY.etc.IMA.imafilename{iko}], '-mat' );

EEG = pop_loadset('filename',IMA.subjfilename{1},'filepath',IMA.subjfilepath{1});


if strcmp(g.stretch_spectra, 'on') && ~isempty(g.peakrange) && isempty(g.targetpeakfreq);
    g.targetpeakfreq = mean(g.peakrange);
elseif strcmp(g.stretch_spectra, 'on') && isempty(g.peakrange);
    fprintf('\n \n cannot stretch spectra without given peak frequency range.... \n')
    fprintf('\n turning off stretch spectra parameter \n \n')
    g.stretch_spectra = 'off';
end



[IMA] = collecttemplates(IMA, EEG, 'peakrange', g.peakrange, 'stretch_spectra', g.stretch_spectra,...
       'targetpeakfreq', g.targetpeakfreq, 'plot_templ', g.plot_templ);


end
end
  