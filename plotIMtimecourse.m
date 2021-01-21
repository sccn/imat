% plot original IC time-frequency decomposition, PC time-frequency backprojection, IM time-frequency backprojection and IM
% weights timecourse
%
% plotIMtimecourse(IMA, vararg)
%
%
% Author: Johanna Wagner, Swartz Center for Computational Neuroscience, UC San Diego, 2019
% adapted from functions written by Julie Onton
%
%
% Example: plot components 1, 5 and 8 with frequency limit 6-40Hz and smoothing 0.1
%           - plot only IC time-frequency decomposition and PC time-frequency backprojections
%
% >>  plotIMtimecourse(IMA, 'comps', [1 5 8],...
%      'frqlim', [6 40], 'plotICtf', 'on', 'plotPCtf', 'on')
%
%
% Example: for all conditions plot components 3, 5 and 6 and IMs 2, 6 and 9 with frequency limit 6-40Hz and smoothing 0.1
%          plot only IM time-frequency backprojections on ICs and IM weight
%           timecourses
%
% >>  pop_plotIMtimecourse_study(IMA, 'comps', [3 5 6], 'factors', [2 6 9],...
%      'frqlim', [6 40], 'plotcond', 'on','smoothing', 0.1, 'plotIMtf', 'on', 'plotIMtime', 'on')
%
%
%
% INPUTS:
% required Inputs:
% IMA - previously saved IMA structure (either by running pop_runima or pop_runima_study)
% EEG - EEG structure of associated EEG dataset
% optional Inputs:
% comps - [vector] independent components to plot - if empty plots all
%                      the independent components for a subject
% factors - [vector] IMs to plot - if empty plots all the IMs for a subject% frqlim - [min max] frequency limits for plotting
% frqlim - [min max] frequency limits for plotting
% plotcond - [string] 'off' or 'on' (default off). If IMA over different
%                      conditions was computed plots a vertical line in between conditions and
%                      plots conditions in different colors for IM weights plotting
% smoothing - [integer] how many timepoints to use for moving average for
%                        smoothing. A higher number means more smoothing. default 40
% plotICtf  -  [string]  {'on' 'off'} plot IC time-frequency decomposition, default 'off'
% plotPCtf  -   [string] {'on' 'off'}  plot PC time-frequency backprojection, default 'off'
% plotIMtf  -    [string] {'on' 'off'}  plot IM time-frequency backprojection, default 'off'
% plotIMtime  -  [string] {'on' 'off'}   plot IM weights timecourse, default 'off'



function plotIMtimecourse(IMA, EEG, varargin)

g = finputcheck(varargin, {'comps'     'integer'   []             []; ...
    'factors'       'integer'    []             []; ...
    'frqlim'        'real'       []             []; ...
    'plotcond'       'string'    {'on' 'off'}                       'off';...
    'smoothing'      'integer'   []             [40];...
    'method'         'string'    {'movmean' 'gaussian' 'lowess' 'rlowess'}  'rlowess';... 
    'plotICtf'       'string'    {'on' 'off'}                       'off';...
    'plotPCtf'       'string'    {'on' 'off'}                       'off';...
    'plotIMtf'       'string'    {'on' 'off'}                       'off';...
    'plotIMtime'     'string'    {'on' 'off'}                       'off';...
    }, 'inputgui');
if isstr(g), error(g); end;


if isempty(g.frqlim)
    g.frqlim = IMA.freqlim;
end

if isempty(g.comps)
    g.comps = IMA.complist;
end

if isempty(g.factors)
    g.factors = 1:IMA.npcs;
end

if strcmp(g.plotcond, 'on') || iscell(IMA.timepntCond);
    dataport = IMA.timepntCond;
else
    dataport = {IMA.timepntCond};
end


%% load EEG dataset associated with IMA for plotting of scalpmaps
%EEG = pop_loadset('filename',IMA.subjfilename{1},'filepath',IMA.subjfilepath{1});

times = IMA.timevec/1000; % transform timevector to seconds
freqvec = IMA.freqvec;
freqscale = IMA.freqscale;
meanpwr = IMA.meanpwr;

sph = IMA.sph;
wts = IMA.wts;
PCact = IMA.pc; % PC spectral backprojection
origspecdat = IMA.timefreq; % original timefrequency matrix containing tf maps of all ICS time x spectra*ICs
speceig = IMA.eigvec; % PC backprojection in time

ws = wts*sph;
winv = pinv(ws);
activations = ws*PCact; % template spectra
specwts = speceig*winv;
winv = specwts; % template timecourse  (overwrite ICA winv with ICA/PCA winv)
clear speceig specwts


%% scale activations
[valAct, indAct] =  max(abs(activations)');
for onj = 1:size(activations,1);
    if activations(onj,indAct(onj))>=0;
        maxval = 1;
    else
        maxval = -1;
    end
    activations(onj,:) = activations(onj,:)*maxval;
    polarity(onj) = maxval;
    winv(:,onj) = winv(:,onj)*maxval;
end




rawbkgd = [.89 .89 .89]; % .89 .89 .89 color of raw data background, [] to turn off
pcabkgd = [.75 .75 .75];%  .75 .75 .75color of pca reduced data background
lnwdth = 2.2;%2.2

cols(1,:) = [0 .8 .8];
cols(2,:) = [1 0 1];
cols(3,:) = [1 0 0];
cols(4,:) = [0 0 1];
cols(5,:) = [0 1 0];
cols(6,:) = [1 .4 0]; % orange
cols(7,:) = [.6 0 .6]; % purple
cols(8,:) = [1 .4 .4]; %light red
cols(9,:) = [.4 .6 0]; % army green
cols(10,:) = [1 .8 0]; % mustard

cols2(1,:) = round([0 .8 .8]*0.9,1);
cols2(2,:) = round([1 0 1]*0.9,1);
cols2(3,:) = round([1 0 0]*0.9,1);
cols2(4,:) = round([0 0 1]*0.9,1);
cols2(5,:) = round([0 1 0]*0.9,1);
cols2(6,:) = round([1 .4 0]*0.9,1); % orange
cols2(7,:) = round([.6 0 .6]*0.9,1); % purple
cols2(8,:) = round([1 .4 .4]*0.9,1); %light red
cols2(9,:) = round([.4 .6 0]*0.9,1); % army green
cols2(10,:) = round([1 .8 0]*0.9,1); % mustard





taglist = {'line1' 'line2' 'line3' 'line4' 'line5' 'line6' 'line7'...
    'line7' 'line8' 'line9' 'line10' 'line11' 'line12' 'line13' 'line14' 'line15'};

icadefs;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% check how many rows/columns are needed for subplots


if length(g.comps) == 2 || length(g.comps) == 3 || length(g.comps) == 4;
    row = length(g.comps);
    col = 3;
else
    row = ceil(sqrt(length(g.comps)*3));
    col = ceil(sqrt(length(g.comps)*3));
    
    rep = 1;
    while rep == 1;
        if col ==3 | col==6 | col ==9 | col==12;
            rep = 0;
        else
            col = col+1; %
        end;
    end;
    rep = 1;
    while rep == 1;
        if length(g.comps)*3 <= (row-1)* col;
            row = row -1;
        else
            rep = 0;
        end;
    end
end;

%row = round(sqrt(length(g.comps))); col = ceil(sqrt(length(g.comps)*3)); % determine how many rows and columns for subplot


cols = hsv(length(g.factors)); % determine colors for IM traces

fr = find(freqvec >= g.frqlim(1) & freqvec <= g.frqlim(2)); % check index for frequency range
freqs = freqvec(fr);
pl = 1;


if strcmp(g.plotICtf, 'on');
    
    % get scale for colorbar
    plotdatascale = [];
    for cp = 1:length(g.comps);
        rcp = find(g.comps(cp) == IMA.complist);
        plotdata= origspecdat(:,length(freqvec)*(rcp-1)+1:length(freqvec)*rcp);
        plotdatascale = [plotdatascale plotdata(:,fr)]; % select frequency range to plot
    end
    % minl = min(plotdata(:));
    maxl = max(abs(plotdatascale(:)))/3; %% changed colorscale min max
    %% plot time frequency map of original spectra
    figure;
    for cp = 1:length(g.comps);
        rcp = find(g.comps(cp) == IMA.complist);
        % plot scalpmap of IC
        sbplot(row,col,[pl]);
        topoplot(EEG.icawinv(:,g.comps(cp)),EEG.chanlocs(1:size(EEG.icawinv,1)),'electrodes','off');
        title(['IC ' int2str(g.comps(cp))]);
        set(gca,'fontsize',14);
        pl = pl+1;
        
        %select tf map of IC from origspecdata and add IC mean spectrum
        %plotdata = origspecdat(:,length(freqvec)*(rcp-1)+1:length(freqvec)*rcp) + repmat(meanpwr(rcp,:),[size(origspecdat(:,:),1) 1]);
        plotdata = origspecdat(:,length(freqvec)*(rcp-1)+1:length(freqvec)*rcp);
        plotdata = plotdata(:,fr); % select frequency range to plot
        % minl = min(plotdata(:));
        %maxl = max(abs(plotdata(:)))/3; %% changed colorscale min max
        sbplot(row,col,[pl pl+1]);
        imagesc(times,freqs,plotdata'); % plot original tf map for IC
        set(gca,'YDir','normal');
        set(gca,'YScale','log');
        set(gca,'ytick',[10 20 40 80 120])
        ylim([freqs(1) freqs(end)])
        caxis([-maxl maxl]);
        j=jet;
        colormap(j);
        colorbar;
        if strcmp(g.plotcond, 'on') % add separating vertical lines between conditions
            for lux = 1:length(dataport);
                line([times(dataport{lux}(end)),times(dataport{lux}(end))],[freqs(1) freqs(end)], 'Color','k', 'LineWidth', lnwdth); hold on
            end
        end
        if pl == (row-1)*col+2
            xlabel('Time (sec)'); ylabel('Frequency (Hz)');
        elseif pl > (row-1)*col+1
            xlabel('Time (sec)');
        end;
        set(gca,'fontsize', 12);
        pl = pl+2;% advance for scalp map and spectra
    end
    ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0  1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
    text(0.4, 0.98,'IC spectogram', 'fontsize', 20);
    if col >4
        set(gcf,'Position',[100 300 1400 900]);
    end
    set(gcf,'PaperOrientation','landscape');  set(gcf,'PaperPosition',[0.25 0.25 10.5 8]);
    set(gcf,'color',BACKCOLOR);
end

%%%%%%%%%%%%%%  Plot PCA backproj %%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(g.plotPCtf, 'on');
    % get scale for colorbar
    plotdatascale = [];
    for cp = 1:length(g.comps);
        rcp = find(g.comps(cp) == IMA.complist);
        pcaproj = winv(:,:)*activations(:,length(freqvec)*(rcp-1)+1:length(freqvec)*rcp); % select IC tf maps from activations
        plotdatascale = [plotdatascale pcaproj(:,fr)]; % select frequency range to plot
    end
    % minl = min(plotdata(:));
    maxl = max(abs(plotdatascale(:)))/2; %% changed colorscale min max
    figure;
    pl = 1;
    
    for cp = 1:length(g.comps);
        rcp = find(g.comps(cp) == IMA.complist);
        
        % plot scalpmap of IC
        sbplot(row,col,[pl]);
        topoplot(EEG.icawinv(:,g.comps(cp)),EEG.chanlocs(1:size(EEG.icawinv,1)),'electrodes','off');
        title(['IC ' int2str(g.comps(cp))]);
        set(gca,'fontsize',14);
        pl = pl+1;
        
        % plot PC backprojection timefrequency map
        pcaproj = winv(:,:)*activations(:,length(freqvec)*(rcp-1)+1:length(freqvec)*rcp); % select IC tf maps from activations
        %pcaproj = pcaproj + repmat(meanpwr(rcp,:),[size(pcaproj,1) 1]); % add mean power of IC
        plotdata = pcaproj(:,fr); % select frequency range to plot
        %     minl = min(plotdata(:));
        %     maxl = max(plotdata(:));
        % maxl = max(abs(plotdata(:)))/2; %% changed colorscale min max
        sbplot(row,col,[pl pl+1]);
        imagesc(times,freqs,plotdata'); % plot tf maps
        set(gca,'YDir','normal');
        set(gca,'YScale','log');
        set(gca,'ytick',[10 20 40 80 120])
        ylim([freqs(1) freqs(end)])
        caxis([-maxl maxl]);
        j=jet;
        colormap(j);
        colorbar;
        if strcmp(g.plotcond, 'on') % add separating vertical lines between conditions
            for lux = 1:length(dataport);
                line([times(dataport{lux}(end)),times(dataport{lux}(end))],[freqs(1) freqs(end)], 'Color','k', 'LineWidth', lnwdth); hold on
            end
        end
        if pl == (row-1)*col+2
            xlabel('Time (sec)'); ylabel('Frequency (Hz)');
        elseif pl > (row-1)*col+1
            xlabel('Time (sec)');
        end;
        set(gca,'fontsize', 12);
        pl = pl+2;% advance for scalp map and spectra
    end
    ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0  1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
    text(0.4, 0.98,'Summed IM backprojection', 'fontsize', 20);
    if col >4
        set(gcf,'Position',[100 300 1400 900]);
    end
    set(gcf,'PaperOrientation','landscape');  set(gcf,'PaperPosition',[0.25 0.25 10.5 8]);
    set(gcf,'color',BACKCOLOR);
    axcopy
end

%% plot IM timefrequency maps
if strcmp(g.plotIMtf, 'on');
    
    for tpp = 1:length(g.factors);
        tp = g.factors(tpp);
        figure;
        pl = 1;
        
        % get scale for colorbar
        plotdatascale = [];
        for cp = 1:length(g.comps);
            rcp = find(g.comps(cp) == IMA.complist);
            backproj = winv(:,tp)*activations(tp,:); % backproj curr IM
            onebkprj = backproj(:,length(freqvec)*(rcp-1)+1:length(freqvec)*rcp);
            plotdatascale = [plotdatascale onebkprj(:,fr)]; % select frequency range to plot
        end
        % minl = min(plotdata(:));
        maxl = max(abs(plotdatascale(:)))/2; %% changed colorscale min max
        
        for cp = 1:length(g.comps);
            rcp = find(g.comps(cp) == IMA.complist);
            
            % plot scalpmap of IC
            sbplot(row,col,[pl]);
            topoplot(EEG.icawinv(:,g.comps(cp)),EEG.chanlocs(1:size(EEG.icawinv,1)),'electrodes','off');
            title(['IC ' int2str(g.comps(cp))]);
            set(gca,'fontsize',14);
            pl = pl+1;
            
            backproj = winv(:,tp)*activations(tp,:); % backproj curr IM
            onebkprj = backproj(:,length(freqvec)*(rcp-1)+1:length(freqvec)*rcp);
            plotdata = onebkprj(:,fr);
            %         minl = min(plotdata(:));
            %         maxl = max(plotdata(:));
            % maxl = max(abs(plotdata(:)))/2; %% changed colorscale min max
            sbplot(row,col,[pl pl+1]);
            imagesc(times,freqs,plotdata');
            set(gca,'YDir','normal');
            set(gca,'YScale','log');
            set(gca,'ytick',[10 20 40 80 120])
            ylim([freqs(1) freqs(end)])
            caxis([-maxl maxl]);
            j=jet;
            colormap(j);
            colorbar;
            if strcmp(g.plotcond, 'on') % add separating vertical lines between conditions
                for lux = 1:length(dataport);
                    line([times(dataport{lux}(end)),times(dataport{lux}(end))],[freqs(1) freqs(end)], 'Color','k', 'LineWidth', lnwdth); hold on
                end
            end
            if pl == (row-1)*col+2
                xlabel('Time (sec)'); ylabel('Frequency (Hz)');
            elseif pl > (row-1)*col+1
                xlabel('Time (sec)');
            end;
            %             if cp == 1;
            %                 title(['IM ' num2str(g.factors(tpp))], 'fontsize', 20);
            %             end
            set(gca,'fontsize', 12);
            
            pl = pl+2;% advance for scalp map and spectra
        end
        ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0  1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
        text(0.4, 0.98,['IM spectral weights backprojection IM ' num2str(g.factors(tpp))], 'fontsize', 16);
        if col >4
            set(gcf,'Position',[100 300 1400 900]);
        end
        set(gcf,'PaperOrientation','landscape');  set(gcf,'PaperPosition',[0.25 0.25 10.5 8]);
        set(gcf,'color',BACKCOLOR);
        axcopy
    end
end

%% plot IM timecourses

%
if strcmp(g.plotIMtime, 'on');
    
    if length(g.factors) == 2 || length(g.factors) == 3 || length(g.factors) == 4;
        row = length(g.factors);
        col = 3;
    else
        row = ceil(sqrt(length(g.factors)*3));
        col = ceil(sqrt(length(g.factors)*3));
        
        rep = 1;
        while rep == 1;
            if col ==3 | col==6 | col ==9 | col==12;
                rep = 0;
            else
                col = col+1; %
            end;
        end;
        rep = 1;
        while rep == 1;
            if length(g.factors)*3 <= (row-1)* col;
                row = row -1;
            else
                rep = 0;
            end;
        end
    end;
    
    plotdatascale = [];
    for tpp = 1:length(g.factors);
        tp = g.factors(tpp);
        backproj = squeeze(winv(:,tp));
        % smooth timecourses using a lowpass butterworth filter
        %         clear aa bb
        %         [bb,aa] = butter (2, 0.1*g.smoothing./(round(times(end)/length(times)*10)/2), 'low');
        %         backproj = filtfilt(bb,aa,(double(backproj)));
        plotdataIM = [];
        for lux = 1:length(dataport);
            plotdataIM = [plotdataIM smoothdata(double(backproj(dataport{lux})),g.method,g.smoothing)']; % ,'movmedian',20
        end
        plotdatascale = [plotdatascale plotdataIM]; % select frequency range to plot
    end
    % minl = min(plotdata(:));
    maxl = max(abs(plotdatascale(:))); %% changed colorscale min max
    
    plotdataIM = [];
    %row = ceil(length(g.factors));
    figure;
    pl = 1;
    for tpp = 1:length(g.factors);
        tp = g.factors(tpp);
        
        % plot histogram of template
        sbplot(row,col,[pl]);
        hist(winv(:,tp),75);hold on;
        title(['IM ' num2str(g.factors(tpp))]);
        set(gca,'fontsize',16);
        plot([0 0],[get(gca,'ylim')],'r-');
        set(gca,'yticklabel',[]);   set(gca,'xticklabel',[]);
        
        pl = pl+1;
        
        backproj = squeeze(winv(:,tp));
        
        % smooth timecourses using a lowpass butterworth filter
        %         clear aa bb
        %         [bb,aa] = butter (2, 0.1*g.smoothing./(round(times(end)/length(times)*10)/2), 'low');
        %         backproj = filtfilt(bb,aa,(double(backproj)));
        %
        % plot IM timecourse
        sbplot(row,col,[pl pl+1]);
        for lux = 1:length(dataport);
            
            plotdataIM = smoothdata(double(backproj(dataport{lux})),g.method ,g.smoothing); % ,'movmedian',20
            line([times(dataport{lux}(1)),times(dataport{end}(end))],[0,0], 'Color','k', 'LineWidth', lnwdth); hold on
            if strcmp(g.plotcond, 'on') % add separating vertical lines between conditions
                line([times(dataport{lux}(end)),times(dataport{lux}(end))],[min(backproj)-1 max(backproj)+1], 'Color','k', 'LineWidth', lnwdth); hold on
            end
            % ph1 = plot(times(dataport{lux}),backproj(dataport{lux}),'LineWidth', lnwdth,'Color',cols(tpp,:));hold on
            ph1 = plot(times(dataport{lux}),plotdataIM,'LineWidth', lnwdth,'Color',cols(tpp,:));hold on
            
            set(ph1,'color',cols(lux,:));
            ph12 = gcf;
            ph12.Children(1).Children(1).Tag = taglist{lux};
            if ~isempty(IMA.condition)
                legendInfo{lux} = ['Cond ' IMA.condition{lux}];
            else
                legendInfo{lux} = ['Cond ' num2str(lux)];
            end
            if pl == (row-1)*col+2;
                xlabel('Time (sec)');
            elseif pl > (row-1)*col+1
                xlabel('Time (sec)');
            end;
        end
        
        xlim([times(dataport{1}(1)),times(dataport{end}(end))]);
        ylim([-maxl maxl]);
        set(gca,'fontsize', 12);
        pl = pl+2;% advance for scalp map and spectra
        if strcmp(g.plotcond, 'on')
            if tpp == 1;
                hb = [];
                for loki = 1:length(IMA.ntrials);
                    
                    hb = [hb findobj(gcf,'tag',taglist{loki})];
                end
                hlegend = legend(hb,legendInfo);
                set(hlegend, 'Location', 'best');
            end
        end
    end
    ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0  1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
    text(0.4, 0.98,'IM weight backprojection', 'fontsize', 20);
    if col >4
        set(gcf,'Position',[100 300 1400 900]);
    end
    icadefs;
    set(gcf,'PaperOrientation','landscape');  set(gcf,'PaperPosition',[0.25 0.25 10.5 8]);
    set(gcf,'color',BACKCOLOR);
    axcopy
end

