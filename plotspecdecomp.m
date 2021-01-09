%
% Plots IMA results: plots IMA templates in a matrix for all IMs and
% ICs separately, or all ICs for each IM superimposed on a single axis, or
% each IC scalp map with super imposed IM templates on a single axis
%
%
% plotspecdecomp(IMA, varargin)
%
% Author: Johanna Wagner, Swartz Center for Computational Neuroscience, UC San Diego, 2019
% adapted from a function written by Julie Onton
%
%
%% Example: plot all ICs and IMs separately for IMA decomposition of a subject
% >> pop_plotspecdecomp(IMA, 'plottype', 'comb')
%
%
% Example: plot all ICs with a selection of superimposed IMs for IMA decomposition of a subject
% >> pop_plotspecdecomp(IMA, EEG, 'plottype', 'ims', 'factors', [1 3 4 6 7], 'maps', 'on')
%
%
%
% INPUTS
% IMA - previously saved IMA structure (created either by running pop_runima or pop_runima_study)
% EEG - EEG structure of associated EEG file
% comps - independent components to plot
% factors - IMs to plot
% frqlim - frequency limits for plotting
% plottype -- ['ims', 'ics' 'comb'] 'ics' will plot a single axis with all comps
%             superimposed.
%             'ims' will show each IC scalp map with super imposed IM templates.
%             'comb' will plot all ICs and IMs plotted separately.
% maps -- ['on','off'] if 'on', then will plot scalp maps.


function plotspecdecomp(IMA, EEG, varargin)

%% check inputs
g = finputcheck(varargin, { 'comps'     'integer'   []             [IMA.complist]; ...
    'factors'       'integer'    []             [1:IMA.npcs]; ...
    'frqlim'        'real'       []             [IMA.freqlim]; ...
    'freqscale'     'string'     {'log' 'linear'}         'log';...
    'plottype'      'string'     {'ics' 'ims' 'comb'}           'comb'; ...
    'maps'          'string'     {'on', 'off'}  'on';...
    }, 'inputgui');
if isstr(g), error(g); end;

plotbackproj = 'off'; % for single IM plotting, plots all trials vs activations
maxrows = 11; % max # of rows to plot before starting new fig
lnwdth = 2;


nlim = []; mlim = [];


%% load EEG dataset associated with IMA for plotting of scalpmaps
% EEG = pop_loadset('filename',IMA.subjfilename{1},'filepath',IMA.subjfilepath{1});

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

    icadefs;

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

scaling = rms(activations,2);
activations = activations./repmat(scaling,1,size(activations,2));

clear wts sph ws icamatall


%% plotting function

fr = find(freqvec >= g.frqlim(1) & freqvec <= g.frqlim(end)); % find index of frequency vector inside freq limit


minl = min(min(activations(g.factors,:)))-abs(min(min(activations(g.factors,:))))*.01;
maxl = max(max(activations(g.factors,:)))+abs(max(max(activations(g.factors,:))))*.01;

%% plot superimposed IC templates for each IM

if strcmp(g.plottype,'ics') % superimpose ic templates
    figure;row = round(sqrt(length(g.factors))); col = ceil(sqrt(length(g.factors))); % determine how many rows and columns for subplot
    cols = lines(length(g.comps)); pl = 1;
    for tpp = 1:length(g.factors)
        tp = g.factors(tpp);
        sbplot(row,col,pl)
        % plot superimposed IC template spectra for each specified IM
        for cp = 1:length(g.comps)
            rcp = find(g.comps(cp) == IMA.complist);
            if strcmp(g.freqscale,'log') % log spaced
                ph = semilogx(freqvec,activations(tp,length(freqvec)*(rcp-1)+1:length(freqvec)*rcp), 'LineWidth', 2);hold on
                set(gca,'FontSize',12)
                set(gca,'xtick',[3 6 10 20 40 80])
                xlim([g.frqlim(1) g.frqlim(end)])
                set(ph,'color',cols(cp,:));
            else % otherwise linear
                ph = plot(freqvec,activations(tp,length(freqvec)*(rcp-1)+1:length(freqvec)*rcp),'linewidth',lnwdth);
                xlim([g.frqlim(1) g.frqlim(end)])
                hold on;
                set(ph,'color',cols(cp,:));
            end;
        end;
        set(gca,'ylim',[minl maxl]); title(['IM ',int2str(g.factors(tpp))]);
        set(gca,'box','off');
        set(gca,'xgrid','on');
        if pl == (row-1)*col+1
            xlabel('Frequency (Hz)'); ylabel('Relative Power');
        elseif pl > (row-1)*col+1
            xlabel('Frequency (Hz)');
        end;
        if pl <= col*(row-1)
            set(gca,'xticklabel',[]);
        end;pl = pl+1;
    end; 
    h90 = textsc(['Superimposed IC Templates for single IMs'],'title');
    set(h90, 'FontSize',20)
    set(gcf,'Position',[100 300 1400 900]);
    set(gcf,'PaperOrientation','landscape');  set(gcf,'PaperPosition',[0.25 0.25 10.5 8]);    
    set(gcf,'color',BACKCOLOR);
    axcopy

    %% plot superimposed IM templates for each IC
    
elseif strcmp(g.plottype,'ims')  % superimpose IM templates for each IC
    figure;
    row = round(sqrt(length(g.comps)*2)); % check how many comlumns and rows for subplot
    col = ceil(sqrt(length(g.comps)*2));
    if mod(col,2) == 1
        col = col+1; row = row-1;
    end;
    cols = lines(length(g.factors));
    pl = 1;
    
    % plot superimposed IC template spectra for each specified IM
    
    for cp = 1:length(g.comps)
        rcp = find(g.comps(cp) == IMA.complist);
        sbplot(row,col,pl); pl = pl+1;
        topoplot(EEG.icawinv(:,g.comps(cp)),EEG.chanlocs(EEG.icachansind),'electrodes','off');
        set(gca,'fontsize',20);  title(['IC ',int2str(g.comps(cp))]);
        sbplot(row,col,pl);
        for tpp = 1:length(g.factors)
            tp = g.factors(tpp);
            if strcmp(g.freqscale,'log') % log spaced
                ph = semilogx(freqvec,activations(tp,length(freqvec)*(rcp-1)+1:length(freqvec)*rcp), 'LineWidth', 2);hold on
                set(gca,'FontSize',12)
                set(gca,'xtick',[3 6 10 20 40 80])
                xlim([g.frqlim(1) g.frqlim(end)])
                set(ph,'color',cols(tpp,:));
            else % otherwise linear
                ph = plot(freqvec,activations(tp,length(freqvec)*(rcp-1)+1:length(freqvec)*rcp),'linewidth',lnwdth);
                xlim([g.frqlim(1) g.frqlim(end)])
                hold on;
                set(ph,'color',cols(tpp,:));
            end;
        end;
        set(gca,'ylim',[minl maxl]); 
        set(gca,'xgrid','on');
        if cp <= round(col/2)
            title(['IM templates'],'fontsize',20);
        end
        if pl == (row-1)*col+2
            xlabel('Frequency (Hz)'); ylabel('Relative Power');
        elseif pl > (row-1)*col+1
            xlabel('Frequency (Hz)');
        end;
        if pl <= col*(row-1)
            set(gca,'xticklabel',[]);
        end;
        pl = pl+1;
    end;
    h90 = textsc(['Superimposed IM templates for single ICs'],'title');
    set(h90, 'FontSize',20)
    set(gcf,'Position',[100 300 1400 900]);
    set(gcf,'PaperOrientation','landscape');  set(gcf,'PaperPosition',[0.25 0.25 10.5 8]);
    set(gcf,'color',BACKCOLOR);
    axcopy

    
    %% plot IM templates separately for each IC
else
    figure;row = length(g.factors)+1;
    if row > maxrows  % check how many columns and rows for subplot
        row = round(row/2);
        if row > maxrows
            row = maxrows;
        end;
    end;
    col = length(g.comps)+1;
    if strcmp(g.maps,'on')% plot scalp maps if requested
        pl = 2;
        for cp = 1:length(g.comps)
            sbplot(row,col,pl)
            topoplot(EEG.icawinv(:,g.comps(cp)),EEG.chanlocs(EEG.icachansind),'electrodes','off','plotrad',.7); pl = pl+1;
            set(gca,'fontsize',14);  title(['IC ' int2str(g.comps(cp))]);
        end;
    else
        pl = 1;
    end;
    
    for tpp = 1:length(g.factors)
        tp = g.factors(tpp);
        if pl == row*col+1
            set(gcf,'Position',[100 300 1400 900]);
            set(gcf,'PaperOrientation','landscape');  set(gcf,'PaperPosition',[0.25 0.25 10.5 8]);
            set(gcf,'color',BACKCOLOR);
            axcopy
            if isempty(EEG.subject) %% plot title
                ph=textsc(['Independent Modulators'],'title');
            else
                ph=textsc(['Independent Modulators'],'title');
            end
            set(ph,'fontsize',20);
            figure;
            if strcmp(g.maps,'on')% plot scalp maps if requested
                pl = 2;
                for cp = 1:length(g.comps)
                    sbplot(row,col,pl)
                    topoplot(EEG.icawinv(:,g.comps(cp)),EEG.chanlocs(EEG.icachansind),'electrodes','off'); pl = pl+1;
                    set(gca,'fontsize',16);  title(['IC' int2str(g.comps(cp))]);
                end;
            else
                pl = 1;
            end;
        end;
        % plot IM histogram
        sbplot(row,col,pl)
        hist(winv(:,tp),75);pl = pl+1;hold on;
        set(gca,'fontsize',7);
        plot([0 0],[get(gca,'ylim')],'r-');
        set(gca,'yticklabel',[]);   set(gca,'xticklabel',[]);
        title(['IM ',int2str(tp)], 'fontsize',12);
        
        % plot template spectra
        for cp = 1:length(g.comps)
            rcp = find(ismember(IMA.complist,g.comps(cp)));
            sbplot(row,col,pl);
            
            if strcmp(g.freqscale,'log')  % log spacing
                ph = semilogx(freqvec,activations(tp,length(freqvec)*(rcp-1)+1:length(freqvec)*rcp), 'LineWidth', 2, 'Color','b');hold on
                set(gca,'FontSize',12)
                set(gca,'xtick',[3 10 30 80])
                xlim([g.frqlim(1) g.frqlim(end)]);
                pl = pl+1;hold on;
            else % otherwise linear
                plot(freqvec,activations(tp,length(freqvec)*(rcp-1)+1:length(freqvec)*rcp),'LineWidth', 2); pl = pl+1;hold on;
                xlim([g.frqlim(1) g.frqlim(end)]);
            end;
            set(gca,'ylim',[minl maxl]);
            set(gca,'ytick',[ceil(minl) 0 floor(maxl)]);
            set(gca,'xgrid','on');
            set(gca,'fontsize',7);set(gca,'box','off');
            set(gca,'ticklength',[.03 .03]);
            plot([get(gca,'xlim')],[0 0],'r-');
            if cp == round(length(g.comps)/2)
                xlabel('Frequency (Hz)')
            end
            set(gca,'fontsize', 12)
            if pl <= (row-1)*col+1
                if tpp ~= length(g.factors)
                    set(gca,'xticklabel',[]);
                    set(gca,'yticklabel',[]);
                    xlabel('')
                elseif tpp ~= length(g.factors) & cp ~= 1;
                    set(gca,'yticklabel',[]);
                end;
            end;
            if ~strcmp(g.maps,'on') & pl <= (col+1)
                title(int2str(g.comps(cp)));
            end;
        end;
    end;
    set(gcf,'Position',[100 300 1400 900]);
    set(gcf,'PaperOrientation','landscape');  set(gcf,'PaperPosition',[0.25 0.25 10.5 8]);
    if isempty(EEG.subject)
        ph=textsc(['Independent Modulators'],'title');
    else
        ph=textsc(['Independent Modulators'],'title');
    end
    set(ph,'fontsize',20);
    set(gcf,'color',BACKCOLOR);
    axcopy

    
end;
