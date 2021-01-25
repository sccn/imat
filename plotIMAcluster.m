% plots spctral templates, dipole density and scalpmaps for IM clusters
%
% plotIMAcluster(clustidx, varargin);
%
% %
% Author: Johanna Wagner, Swartz Center for Computational Neuroscience, UC San Diego, 2019
% adapted from a function written by Julie Onton
%
% Example plotting cluster numbers 1 and 2 and plotting spectra with log
% frequency scale
% Hz and log frequency scale
%  >>  plotIMAcluster(clustidx, 'plottemplates', 'on', 'templates', templates,...
%                      'freqvec', freqvec, 'clust', [1 2], 'freqscale', 'log')
%
% Example plotting subclusters and scalpmaps
%  >>  plotIMAcluster(clustidx, 'plotscalpmaps', 'on', 'scalpmaps', scalpmaps)
%
% Example plotting subclusters, dipoles and spectral templates
%  >>  plotIMAcluster(clustidx, 'plotsubclusters', 'on', 'plottemplates', 'on',...
%                     'templates', templates, 'freqvec', freqvec, 'plotscalpmaps', 'on', 'scalpmaps', scalpmaps)
%
%
% INPUTS:
% required Inputs:
% clustidx - matrix of indices of subjects and ICs in each cluster - output of
%             clusterIMAtemplates, also saved in STUDY.etc.IMA.clustidx
%             clustix contains the indexes of [subjectID, IM, IC, cluster, subcluster]
%
% optional Inputs:
% plotclust - 'string' 'on' or 'off' - default 'on' - may be turned off
%               when wanting to plot only subclusters
% clust - [vector] 'vector of integers' vector of cluster indices to plot
%                  if empty plot all clusters
% freqscale -  'string' frequency scale to plot spectra 'linear' or 'log'
%                       default 'log'
% plotsubclusters - 'string' 'on' or 'off' specifies whether to plot
%                    subclusters - plots only the subclusters previously
%                    defined in pop_subclusterIMAtemplates - default 'off'
% plotscalpmaps - 'string' 'on' or 'off' specifies whether to plot
%                  scalpmaps - default 'off'
% templates - matrix of spectral templates (ICs x spectra) from all subjects and ICs
%              that were included in the clustering (see pop_plotIMAcluster on how the
%              structure is constructed from single subject files) - is required if 'plottemplates is 'on'
% dipsources - structure of dipole locations - as saved in collecttemplates
%              and pop_collecttemplates for single subjects - note that here dipsources
%              should be a structure of dipole locations from all subjects and ICs
%              that were included in the clustering (see pop_plotIMAcluster on how the
%              structure is constructed from single subject files) - is required if 'plotdipsources is 'on'
% scalpmaps - matrix of scalpmaps from all subjects and ICs that were included in the clustering (see pop_plotIMAcluster on how the
%              structure is constructed from single subject files)- is required if 'plotscalpmaps is 'on'
% freqvec - frequency vector of spectra - this is the original frequency
%           vector saved in IMA.freqvec, or adapted if more restrictive frequency
%           limits are chosen in pop_plotIMAcluster - required if plotting
%           spectral templates

%
%
% OUTPUTS: plots of dipoledensity, spectra and scalpmaps of clusters and subclusters
%

function plotIMAcluster(clustidx, varargin);


g = finputcheck(varargin, {'clust'      'integer'     []         [];...
    'plotclust'      'string'   {'on' 'off'}       'on';...
    'freqscale'        'string'    {'linear' 'log'}       'log';...
    'freqvec'        'integer'    []      [];...
    'plotsubclusters'     'string'     {'on' 'off'}         'off';...
    'plottemplates'     'string'    {'on' 'off'}      'off';...
    'plotdipsources'     'string'    {'on' 'off'}      'off';...
    'plotscalpmaps'   'string'     {'on' 'off'}         'off';...
    'templates'      'real'         []          [];...
    'dipsources'     'struct'    {}          {};...
    'scalpmaps'     'real'     []    [];...
    }, 'inputgui');
if isstr(g), error(g); end;


if isempty(g.clust);
    indnonzero = find(clustidx(:,4));
    g.clust = unique(clustidx(indnonzero,4));
end

denswt = zeros(1,0);

row = round(sqrt(length(g.clust))); col = ceil(sqrt(length(g.clust)));
% col = round(length(g.clust)/2);
% if length(g.clust) == 1;
%     col  = 1;
% elseif round(length(g.clust)/2) == 1;
%     col  = col+1;
% end

if strcmp(g.plotclust, 'on');
    if strcmp(g.plotdipsources, 'on') && ~isempty(g.dipsources);
        
        for dind = 1:length(g.clust);
            
            indclus = find(clustidx(:,4) == g.clust(dind));
            dipoles_clus = g.dipsources(indclus);
            figure;
            optdipplot = {dipoles_clus,'gui','off','image','mri','coordformat','spherical','normlen','on'};
            [dens3d mri] = dipoledensity( optdipplot, 'method','alldistance','methodparam',15,'weight',denswt);
            mri3dplot(dens3d,mri); %, 'cmax', 0.08, 'cmin', 0
            tmp = gcf;
            tmp.Name= ['cluster ' num2str(g.clust(dind)) ' (SJs ' num2str(length(unique(clustidx(indclus,1))))...
                '  STs ' num2str(length(clustidx(indclus,2))) ')'];
            t = title(['cluster ' num2str(g.clust(dind)) ' (SJs ' num2str(length(unique(clustidx(indclus,1))))...
                '  STs ' num2str(length(clustidx(indclus,2))) ')'], 'FontSize', 16)
            
        end
    end
    
    if strcmp(g.plotscalpmaps, 'on');
        figure;
        for dind = 1:length(g.clust);
            
            indclus = find(clustidx(:,4) == g.clust(dind));
            scalpmaps_clus = g.scalpmaps(indclus,:,:);
            AVscalp = squeeze(mean(scalpmaps_clus,1));
            
            subplot(row, col,dind);
            toporeplot(AVscalp, 'plotrad', 0.7, 'intrad', 0.5, 'colormap', 'jet');

            t = title(['cluster ' num2str(g.clust(dind)) ' (SJs ' num2str(length(unique(clustidx(indclus,1))))...
                '  STs ' num2str(length(clustidx(indclus,2))) ')'], 'FontSize', 16)
        end
    end
    
    if strcmp(g.plottemplates, 'on') && ~isempty(g.templates)
        
        plotdatascale = [];
        for dind = 1:length(g.clust);
            indclus = find(clustidx(:,4) == g.clust(dind));
            templates_clus = g.templates(indclus,:);
            plotdatascale = [plotdatascale; templates_clus];
        end
        
        maxl = max((plotdatascale(:)))+1; %% changed colorscale min max
        minl = min((plotdatascale(:)))-1;
        
        figure;
        for dind = 1:length(g.clust);
            
            indclus = find(clustidx(:,4) == g.clust(dind));
            templates_clus = g.templates(indclus,:);
            scalpmaps_clus = g.scalpmaps(indclus,:,:);
            AVscalp = squeeze(mean(scalpmaps_clus,1));
            
            subplot(row, col,dind);
            if strcmp(g.freqscale,'log');
                xlog = logspace(log10(g.freqvec(1)), log10(g.freqvec(end)), 8);
                semilogx(g.freqvec,templates_clus', 'LineWidth', 2,'Color','m');hold on
                semilogx(g.freqvec,mean(templates_clus,1)', 'LineWidth', 2,'Color','b');
                set(gca,'FontSize',12);
                set(gca,'xtick',[10 20 40 80]);
                xlim([g.freqvec(1) g.freqvec(end)]);
                % if dind == length(g.clust)
                %  xlabel('Frequency Hz');
                %end
            else
                ph = plot(g.freqvec,templates_clus','m-','linewidth',2); hold on;
                ph = plot(g.freqvec,mean(templates_clus,1)','b-','linewidth',2);
                set(gca,'xlim',[g.freqvec(fr(1)) g.freqvec(end)]);
                realx = get(gca,'xtick'); labelx = get(gca,'xticklabel');
                xlabel('Frequency Hz');
            end;
            set(gca,'ylim',[minl maxl])
            set(gca,'ticklength',[.05 .05]);
            plot([get(gca,'xlim')],[0 0],'k-'); hold on;
            plot([10 10],[get(gca,'ylim')],'g-', 'LineWidth',2); hold on;
            plot([20 20],[get(gca,'ylim')],'g-', 'LineWidth',2); hold on;
            if dind == (row-1)*col+1
                xlabel('Frequency (Hz)'); ylabel('Relative Power');
            elseif dind > (row-1)*col+1
                xlabel('Frequency (Hz)');
            end;
            if dind <= col*(row-1)
                set(gca,'xticklabel',[]);
            end;
            t = title(['cluster ' num2str(g.clust(dind)) ' (SJs ' num2str(length(unique(clustidx(indclus,1))))...
                '  STs ' num2str(length(clustidx(indclus,2))) ')'], 'FontSize', 16)
            
        end
    end
end



if strcmp(g.plotsubclusters, 'on');
    
indnonzero = find(clustidx(:,5));
numsubclus = unique(clustidx(indnonzero,5));
assign_subclus = unique(clustidx(indnonzero,4));

col = round(length(numsubclus)/2);

if length(numsubclus) == 1;
    col  = 1;
elseif round(length(numsubclus)/2) == 1;
    col  = col+1;
end   
    
    
    if strcmp(g.plotdipsources, 'on') && ~isempty(g.dipsources);
        
        
        
            
               indsubclus = find(clustidx(:,5) > 0);
               indclus = unique(clustidx(indsubclus,4));
        for dind = indclus';
            indclusplot = find(clustidx(:,4) == dind);
            subclusplot = unique(clustidx(indclusplot,5));
            
            for dins = subclusplot';
                
                indsubclusplot = find(clustidx(indclusplot,5) == dins);
                
                dipoles_clus = g.dipsources(indsubclusplot);
                figure;
                optdipplot = {dipoles_clus,'gui','off','image','mri','coordformat','spherical','normlen','on'};
                [dens3d mri] = dipoledensity( optdipplot, 'method','alldistance','methodparam',15,'weight',denswt);
                mri3dplot(dens3d,mri); %, 'cmax', 0.08, 'cmin', 0
                tmp = gcf;
                tmp.Name= ['cluster ' num2str(dind) ' subcluster ' num2str(dins)];
            end
        end
    end
    if strcmp(g.plotscalpmaps, 'on') && ~isempty(g.scalpmaps);
        
        indsubclus = find(clustidx(:,5) > 0);
               indclus = unique(clustidx(indsubclus,4));
        for dind = indclus';
            indclusplot = find(clustidx(:,4) == dind);
            subclusplot = unique(clustidx(indclusplot,5));
            figure;
            for dins = subclusplot';
                
                indsubclusplot = find(clustidx(indclusplot,5) == dins);
                
                scalpmaps_clus = g.scalpmaps(indsubclusplot,:,:);
                AVscalp = squeeze(mean(scalpmaps_clus,1));
                subplot(round(length(subclusplot')/2), col,dins);
                toporeplot(AVscalp, 'plotrad', 0.7, 'intrad', 0.5, 'colormap', 'jet');
                title(['cluster ' num2str(dind) ' subcluster ' num2str(dins)], 'FontSize',16);
                
            end
        end
    end
    
    
    if strcmp(g.plottemplates, 'on') && ~isempty(g.templates);
        
        indsubclus = find(clustidx(:,5) > 0);
               indclus = unique(clustidx(indsubclus,4));
        for dind = indclus';
            indclusplot = find(clustidx(:,4) == dind);
            subclusplot = unique(clustidx(indclusplot,5));
            figure;
            for dins = subclusplot';
                
                indsubclusplot = find(clustidx(indclusplot,5) == dins);
                
                templates_clus = g.templates(indsubclusplot,:);
                subplot(round(length(subclusplot')/2), col,dins);
                 
                if strcmp(g.freqscale,'log');
                    %xlog = logspace(log10(g.freqvec(1)), log10(g.freqvec(end)), 8);
                    semilogx(g.freqvec,templates_clus', 'LineWidth', 2,'Color','m');hold on
                    semilogx(g.freqvec,mean(templates_clus,1)', 'LineWidth', 2,'Color','b');
                    set(gca,'FontSize',12);
                    set(gca,'xtick',[10 20 40]);
                    xlim([g.freqvec(1) g.freqvec(end)]);
                    %if dins == length(subclus)
                        xlabel('Frequency Hz');
                    %end
                else
                    ph = plot(g.freqvec,templates_clus','m-','linewidth',2); hold on;
                    ph = plot(g.freqvec,mean(templates_clus,1)','b-','linewidth',2);
                    set(gca,'xlim',[g.freqvec(fr(1)) g.freqvec(end)]);
                   % realx = get(gca,'xtick'); labelx = get(gca,'xticklabel');
                    xlabel('Frequency Hz');
                end;
                set(gca,'ticklength',[.05 .05]);
                plot([get(gca,'xlim')],[0 0],'k-'); hold on;
                plot([10 10],[get(gca,'ylim')],'g-', 'LineWidth',2); hold on;
                plot([20 20],[get(gca,'ylim')],'g-', 'LineWidth',2); hold on;
                title(['cluster ' num2str(dind) ' subcluster ' num2str(dins)], 'FontSize',16);
                
            end
        end
    end
end
