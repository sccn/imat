% 
% stretches single trial spectra to a target frequency using a reference peak
% frenquency i.e. stretches a spectrum to the target frequency 10 using the
% subject specific mean frequency 11.5 calculated from the median template
% peak frequency
% 
% [warpmat, TWseries] = warp_templates(tfdata, evLat, newLat, timevec);
%
% Author: Johanna Wagner, Swartz Center for Computational Neuroscience, UC San Diego, 2019
% adapted from a function written by Julie Onton
%
% Example
% >> [warpmat, Freqtw] = warp_templates(spectemplates, [11.2 11.2 11.2 11.2 11.2], 10, freqvec);
% 
%
% INPUTS
%
% tfdata - [matrix] single trial spectra (trials x spectra)
% evLat - [real] are all the peak frequencies to warp in format (trials x peak)
% newLat - [real] are the new target peak frequencies to which to warp (row vector) 
% timevector - [vector] frequency vector
%
% OUTPUTS:
% TWseries -  the warped templates
% warpmat - the warping matrix with which the original spectra were
%           multiplied to warp the data

function [warpmat, TWseries] = warp_templates(tfdata, evLat, newLat, timevec);

timevecabsS = (abs(timevec(1))+timevec);


clear indexAtMinL;
clear indexAtMinNew;

for k = 1: size(evLat,1);
for i = 1:size(evLat,2);
[minDifferenceValue, indexAtMinL(k,i)] = min(abs(timevecabsS - (abs(timevec(1))+evLat(k,i))));
end
end

for i = 1:size(evLat,2);
[minDifferenceValue, indexAtMinNew(i)] = min(abs(timevecabsS - (abs(timevec(1))+newLat(i))));
end


TWseries = zeros(size(tfdata));


for i = 1:size(evLat,1);
warpmat(i,:,:) = timewarp([1 indexAtMinL(i,:) size(timevec,2)], [1 indexAtMinNew size(timevec,2)]);
TWseries(i,:) = transpose(squeeze(warpmat(i,:,:))*squeeze(tfdata(i,:))');
display(['timewarping...trial...' num2str(i) ' of ' num2str(size(evLat,1))])
end


end