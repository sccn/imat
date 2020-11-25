% takes spectra and a freq vector and stretches the median template peak of all spectral templates to a given frenquency
%
%
% [warped_templs, AV_tempFreq, warpedfreq] = StretchFreqs(spectemplates,freqvec,peaklims,targetpeakfreq);
%
% Author: Johanna Wagner, Swartz Center for Computational Neuroscience, UC San Diego, 2019
% adapted from a function written by Julie Onton
%
%  [warped_templs, AV_tempFreq, warpedfreq] = StretchFreqs(templates,freqvec,[8 14],10);
 
%
% INPUTS:
% spectemplates -- [matrix templates x freqs] - each row is a spectral template - spectral templates to warp
% freqvec--  [vector] frequency vector
% peaklims -- [minfreq maxfreq] min and max for interval in which peak should be 
% targetpeakfreq -- [integer] frequency to which template spectra should be warped to
%
% OUTPUTS:
% warped_templs -- warped templates.
% AV_tempFreq -- median template peak frequency used for warping
% warpedfreq -- median template peak frequency after warping

function [warped_templs, AV_tempFreq, warpedfreq] = StretchFreqs(spectemplates,freqvec,peaklims,targetpeakfreq);
 


%% get median frequency peak from templates for warping

[val, peakfreqind1] = min(abs(freqvec-(peaklims(1)-2)));
[val, peakfreqind2] = min(abs(freqvec-(peaklims(2)+2)));

freqvecpeaklim = freqvec(peakfreqind1:peakfreqind2);
PeakFreqTemplates = [];
for ki= 1:size(spectemplates,1);
[PKS, LOCS] = findpeaks(spectemplates(ki,peakfreqind1:peakfreqind2));
[maxPKS, ind] = max(PKS);
indPK = find(PKS == maxPKS);
PeakFreqTemplates(ki) = freqvecpeaklim(LOCS(indPK));
end

AV_tempFreq = median(PeakFreqTemplates);

%% warp templates to target frequency

 [warpmat, Freqtw] = warp_templates(spectemplates, repmat(AV_tempFreq,size(spectemplates,1),1), targetpeakfreq, freqvec);
       

%% get warped template peak frequencies     

PeakFreqTemplatesWARPED = [];
for ki= 1:size(Freqtw,1);
[PKS, LOCS] = findpeaks(Freqtw(ki,peakfreqind1:peakfreqind2));
[maxPKS, ind] = max(PKS);
indPK = find(PKS == maxPKS);
PeakFreqTemplatesWARPED(ki) = freqvecpeaklim(LOCS(indPK));
end

warpedfreq = median(PeakFreqTemplatesWARPED);

    fprintf('done.\n');
    
 warped_templs = Freqtw;
end

    

