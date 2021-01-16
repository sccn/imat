
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
% >> [IMA] = pop_collecttemplates(STUDY, 'peakfreq', [8 14],'stretch_spectra', 'on', 'targetpeakfreq', 10, 'plot_template', 'on');
%
% Example: chooses all templates whose spectral amplitude are larger than 
% 1/3 of the maximum template amplitude
% >> [IMA] = pop_collecttemplates(STUDY);
%
%
%
% INPUTS:
% required Inputs:
% IMA - previously saved IMA structure (either by running pop_runima or pop_runima_study)
% EEG - EEG structure of associated EEG dataset
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


function [IMA] = collecttemplates(IMA, EEG, varargin);

g = finputcheck(varargin, {'peakrange'        'integer'       []             []; ...
    'stretch_spectra'     'string'     {'on' 'off'}         'off';...
    'targetpeakfreq'      'integer'       []             []; ...
    'plot_templ'             'string'    {'on' 'off'}    'off';...
    }, 'inputgui');
if isstr(g), error(g); end;


if strcmp(g.stretch_spectra, 'on') && ~isempty(g.peakrange) && isempty(g.targetpeakfreq);
    g.targetpeakfreq = mean(g.peakrange);
elseif strcmp(g.stretch_spectra, 'on') && isempty(g.peakrange);
    fprintf('\n \n cannot stretch spectra without given peak frequency range.... \n')
    fprintf('\n turning off stretch spectra parameter \n \n')
    g.stretch_spectra = 'off';
end


%EEG = pop_loadset('filename',IMA.subjfilename{1},'filepath',IMA.subjfilepath{1});

times = IMA.timevec/1000; % transform timevector to seconds
freqvec = IMA.freqvec;
freqscale = IMA.freqscale;
meanpwr = IMA.meanpwr;

sph = IMA.sph;
wts = IMA.wts;
PCact = IMA.pc; % PC spectral backprojection % former ICAmatall
origspecdat = IMA.timefreq; % original timefrequency matrix containing tf maps of all ICS time x spectra*ICs
speceig = IMA.eigvec; % PC backprojection in time

ws = wts*sph;
winv = pinv(ws);
activations = ws*PCact; % template spectra
specwts = speceig*winv;
winv = specwts; % template timecourse  (overwrite ICA winv with ICA/PCA winv)
clear speceig specwts



%% fix ploarity of IMs in the clus_frqlim frequency range
[valActMAX, indAct] =  max(abs(activations)');
for onj = 1:size(activations,1);
    if activations(onj,indAct(onj))>=0;
        polarity(onj,:) =  1;
    else
        polarity(onj,:) = -1;
    end
end

templatesPol = activations.*repmat(polarity,1,size(activations,2));

%% scale templates by RMS
pu = rms(templatesPol');

templates_scaled = templatesPol./repmat(pu',1,size(templatesPol,2));


templatesIMIC = reshape(templates_scaled(:,:)', size(templates_scaled,2)/length(IMA.complist) ,size(templates_scaled,1)*length(IMA.complist))';



IndIMIC = [repelem(1:size(templates_scaled,1),length(IMA.complist))', repmat(IMA.complist,1,size(templates_scaled,1))'];

% select templates in specific frequency range

if ~isempty(g.peakrange)
tempvec = ones(1,length(freqvec))*-1;
%
templ_frqlim = g.peakrange;
[val, freqind1] = min(abs(freqvec-templ_frqlim(1)));
[val, freqind2] = min(abs(freqvec-templ_frqlim(2)));

%window = barthannwin(length(freqvec(freqind1:freqind2)));
window = tukeywin(length(freqvec(freqind1:freqind2)),0.9);
tempvec(freqind1:freqind2) = window;




for inj = 1:size(templatesIMIC,1)
    r=corrcoef(tempvec, templatesIMIC(inj,:));
    correl(inj) = r(1,2);
end
indICs = find(correl>0.4);

else

indICs = 1:size(IndIMIC,1);
    
end

ValMax = max(templatesIMIC(indICs,:)');
[ValMax2, inMax] = max(ValMax);
indKeepIMs = find(ValMax > ValMax2/3);

TempKeepInd = IndIMIC(indICs(indKeepIMs),:);



if strcmp(g.stretch_spectra, 'on');
    
    [warped_templs, AV_tempFreq, warpedfreq] = StretchFreqs(templatesIMIC(indICs(indKeepIMs),:),freqvec,g.peakrange,g.targetpeakfreq);
    TempKeep = warped_templs;
    
else
    
   
    TempKeep = templatesIMIC(indICs(indKeepIMs),:);
end


% collect dipole locations for clustering

dipsources = [];
for temp_idx = 1:size(TempKeepInd,1)
    figure;
    dipsources(1,temp_idx).posxyz = EEG.dipfit.model(TempKeepInd(temp_idx,2)).posxyz;%% index should be IC index
    dipsources(1,temp_idx).momxyz = EEG.dipfit.model(TempKeepInd(temp_idx,2)).momxyz;
    dipsources(1,temp_idx).rv = EEG.dipfit.model(TempKeepInd(temp_idx,2)).rv;
    [h grid_or_val plotrad_or_grid, xmesh, ymesh] = topoplot(EEG.icawinv(:,TempKeepInd(temp_idx,2)),EEG.chanlocs(1:size(EEG.icawinv,1)),'electrodes','off', 'noplot', 'on');
    scalpmaps(temp_idx,:,:) = grid_or_val;
end;
%

if strcmp(g.plot_templ, 'on')
if ~isempty(g.peakrange)
figure;
subplot(2,2,1)
plot(templatesIMIC')
semilogx(freqvec,templatesIMIC', 'LineWidth', 1);hold on
set(gca,'FontSize',12)
set(gca,'xtick',[10 20 40 80 120])
xlim([freqvec(1) freqvec(end)])
title(['All IM templates ' IMA.subj{:}])
subplot(2,2,2)
semilogx(freqvec,tempvec, 'LineWidth', 2,'Color','m');hold on
set(gca,'FontSize',12)
set(gca,'xtick',[10 20 40 80 120])
xlim([freqvec(1) freqvec(end)])
title('Freq range selected')
subplot(2,2,3)
semilogx(freqvec,templatesIMIC(indICs(indKeepIMs),:)', 'LineWidth', 2,'Color','m');hold on
set(gca,'FontSize',12)
set(gca,'xtick',[10 20 40 80 120])
xlim([freqvec(1) freqvec(end)])
title(['Selected templates ' IMA.subj{:}])
xlabel('Frequency (Hz)')
else
figure;
subplot(1,2,1)
plot(templatesIMIC')
semilogx(freqvec,templatesIMIC', 'LineWidth', 1);hold on
set(gca,'FontSize',12)
set(gca,'xtick',[10 20 40 80 120])
xlim([freqvec(1) freqvec(end)])
title(['All IM templates ' IMA.subj{:}])
subplot(1,2,2)
semilogx(freqvec,templatesIMIC(indICs(indKeepIMs),:)', 'LineWidth', 2,'Color','m');hold on
set(gca,'FontSize',12)
set(gca,'xtick',[10 20 40 80 120])
xlim([freqvec(1) freqvec(end)])
title(['Selected templates ' IMA.subj{:}])    
xlabel('Frequency (Hz)')
end
end




IMA.precluster.templates = TempKeep;
IMA.precluster.IMICindex = TempKeepInd;
IMA.precluster.dipsources = dipsources;
IMA.precluster.scalpmaps = scalpmaps;

%% save IMA results to file 
save([IMA.subjfilepath{1} '/' IMA.filename],'IMA');











