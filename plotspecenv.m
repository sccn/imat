% plot envelopes of IC spectra, PC backprojection and IM templates on the same axis
%
% plotspecenv(IMA,varargin);
%
%
% Author: Johanna Wagner, Swartz Center for Computational Neuroscience, UC San Diego, 2019
% adapted from a function written by Julie Onton
%
%
% Example plotting the spectra for a subject for ICs 1 3 5 6 and IMs 1 2 3
% 5 7 with frequency limits 3 to 40 for condition 2
%
% >> plotspecenv(IMA, 'comps', [1 3 5 6], 'factors', [1 2 3 5 7], 'frqlim', [3 40],...
%                 'dataport', IMA.timepntCond{2}, 'meanspec', IMA.meanpwrCond(2,:,:));
%
%
% Example plotting the overall spectra for a subject for IC 1 and IMs 1 and 2
%  with frequency limits 3 to 60, plotting the upper 95th percentile over
%  two conditions or over one condition (if IMA has been computed only over one condition)
%
% >> plotspecenv(IMA, 'comps', [1], 'factors', [1 2], 'frqlim', [3 60], 'plotenv', 'upper');
%
%
% INPUTS:
% IMA - previously saved IMA structure
% EEG - EEG structure of associated EEG file
% comps - [vector] independent components to plot - if empty plots all
%                      the independent components for a subject
% factors - [vector] IMs to plot - if empty plots all the IMs for a subject% frqlim - [min max] frequency limits for plotting
% plotenv - [string] possible values: 'env', 'upper', 'lower'. Either plot envelope of IMs, plot upper 95th percentile,
%           plot lower 95th percentile. (default 'env')
% plotperc - [number between 0 and 1] define a percentile of trials to use to exclude outliers...
%            (default 95th percentile - 0.95)
% setmimax - [min max] set minimum / maximum of spectrum yscale - if empty
%                      uses matlab default settings
% dataport - vector of timepoints to plot - can be either IMA.timepntCond{ib} if IMA has been computed
%            over several conditions using the EEGLAB STUDY structure and conditions are plotted separately
%            or 1:(IMA.ntrials*IMA.ntw_trials) if the overall spectra over conditions is to be
%            plotted, or IMA was only computed over a single condition.
%            default: [1:(IMA.ntrials*IMA.ntw_trials)]
% meanspec - meanspectra of ICs to plot can be either IMA.meanpwrCond(i,:,:) if IMA has been computed
%            over several conditions using the EEGLAB STUDY structure and conditions are plotted separately (see also dataport)
%            or IMA.meanpwr if the overall spectra over conditions is to be
%            plotted, or IMA was only computed over a single condition.
%            default [IMA.meanpwr]
%


function plotspecenv(IMA, EEG, varargin);

g = finputcheck(varargin, { 'comps'     'integer'   []             []; ...
    'factors'          'integer'    []             []; ...
    'frqlim'           'real'       []             []; ...
    'plotenv'          'string'     {'env' 'upper' 'lower'}           'env'; ...
    'plotperc'         'real'       []           [0.95]; ...
    'setminmax'        'real'      []             []; ...
    'dataport'         'real'      []             []; ...
    'meanspec'         'real'      []             []; ...
    }, 'inputgui');
if isstr(g), error(g); end;


if isempty(g.dataport)
    g.dataport = 1:(IMA.ntrials*IMA.ntw_trials);
end

if isempty(g.meanspec)
    g.meanspec = IMA.meanpwr;
end

if isempty(g.frqlim)
    g.frqlim = IMA.freqlim;
end

if isempty(g.comps)
    g.comps = IMA.complist;
end

if isempty(g.factors)
    g.factors = 1:IMA.npcs;
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



%% define colors for plotting
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

taglist = {'line1' 'line2' 'line3' 'line4' 'line5' 'line6' 'line7'...
    'line7' 'line8' 'line9' 'line10' 'line11' 'line12' 'line13' 'line14' 'line15'...
    'line16' 'line17' 'line18' 'line19' 'line20' 'line21' 'line22' 'line23' 'line24'...
    'line25' 'line26' 'line27' 'line28'};


%% check how many columns/rows are needed for subplots
%hfig = figure('Units', 'normalized', 'Position', [0.1149 0.1000 0.7393 0.7790]);
figure;
if length(g.comps) == 2 || length(g.comps) == 3 || length(g.comps) == 4;
    row = length(g.comps);
    col = 3;
else
    row = ceil(sqrt(length(g.comps)*3));
    col = ceil(sqrt(length(g.comps)*3));
    
    rep = 1;
    while rep == 1;
        if col ==3 | col==6 | col==9 | col==12 | col ==18 | col == 24 | col == 30;
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

if length(g.comps) == 1;
    row = 1;
    col = 4;
end




cols = hsv(length(g.factors)); % determine colors for IM traces

% find index for frequency vector limits
fr = find(freqvec >= g.frqlim(1) & freqvec <= g.frqlim(2));
freqs = freqvec(fr);
pl = 1;


%% plotting envelopes of spectras

for cp = 1:length(g.comps)
    rcp = find(g.comps(cp) == IMA.complist);
    
    % plot scalpmaps
    sbplot(row,col,[pl pl]);
    topoplot(EEG.icawinv(:,g.comps(cp)),EEG.chanlocs(1:size(EEG.icawinv,1)),'electrodes','off');
    title(['IC ' int2str(g.comps(cp))]);
    set(gca,'fontsize',14);
    pl = pl+1;
    
    % get spectrum for current IC and add mean IC spectrum
    plotdata = origspecdat(:,length(freqvec)*(rcp-1)+1:length(freqvec)*rcp) + repmat(g.meanspec(rcp,:),[size(origspecdat(:,:),1) 1]);
    plotdata = plotdata(:,fr);
    sbplot(row,col,[pl+0.2 pl+1]);
    
    
    %%%%%%%%%%%%%%  Plot full data %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % remove outlier data by removing spectra of timewindows that exceed predefined percentage (i.e. 0.95 % of the data)
    plotdata2 = zeros(2,length(freqs));
    [tpsort sortidx] = sort(plotdata,1);% sort in timewindow dimension
    if g.plotperc < 1 % shave the outermost timwindows
        cutends = round(size(tpsort,1)-(size(tpsort,1)*g.plotperc));
        tpsort(1:cutends,:) = [];
        tpsort(end-cutends:end,:) = [];
    end;
    mdn = median(tpsort);
    plotdata2 = [tpsort(end,:);tpsort(1,:)];
    
    
    %% plot enevlope of trimmed data
    if ~isempty(rawbkgd)
        if strcmp(freqscale,'linear') %% -- plot as linear spaced frequencies
            f1 = [freqs,[freqs(length(freqs):-1:1)]];
            f2 = [plotdata2(1,:),[plotdata2(2,end:-1:1)]];
            ph = fill(f1,f2,rawbkgd);hold on;
        elseif strcmp(freqscale,'log') %%  --- log-spaced freqs
            f1 = [freqs,[freqs(length(freqs):-1:1)]];
            f2 = [plotdata2(1,:),[plotdata2(2,end:-1:1)]];
            ph = fill(f1,f2,rawbkgd);hold on;
            set(gca,'XScale','log')
            set(gca,'FontSize',12)
            set(gca,'xtick',[10 20 40 80 120])
            xlim([freqs(1) freqs(end)])
            if cp < length(g.comps)
                set(gca,'xticklabel',[]);
            end
        end;
        set(ph,'edgecolor',rawbkgd); %set(ph,'edgealpha',0);
        minl = min(plotdata2(:));
        maxl = max(plotdata2(:));
        clear plotdata plotdata2 envdata
    else
        minl = 0; maxl = 0;
    end;
    
    
    %%%  Plot PCA-reduced data back-projection  %%%%%%%%%%%%%%%%%%%%
    pcaproj = winv(:,:)*activations(:,length(freqvec)*(rcp-1)+1:length(freqvec)*rcp); % get backprojection of current IC
    pcaproj = pcaproj + repmat(g.meanspec(rcp,:),[size(pcaproj,1) 1]); % add mean spectrum
    pcaproj=pcaproj(:,fr); % select frequency range
    
    % remove outlier data by removing spectra of timewindows that exceed predefined percentage (i.e. 0.95 % of the data)
    pcaproj2 = zeros(2,length(freqs));
    [tpsort sortidx] = sort(pcaproj,1);% sort in pc dimension
    if g.plotperc < 1 % shave the outermost timewindows
        cutends = round(size(tpsort,1)-(size(tpsort,1)*g.plotperc));
        tpsort(1:cutends,:) = [];
        tpsort(end-cutends:end,:) = [];
    end;
    mdn = median(tpsort);
    pcaproj2 = [tpsort(end,:);tpsort(1,:)];
    
    
    %% plot the trimmed data
    if ~isempty(pcabkgd) % only plot if given a color
        if strcmp(freqscale,'linear') %% -- plot as linear spaced frequencies
            f1 = [freqs,[freqs(length(freqs):-1:1)]];
            f2 = [pcaproj2(1,:),[pcaproj2(2,end:-1:1)]];
            ph = fill(f1,f2,pcabkgd);hold on;
        elseif strcmp(freqscale,'log') %% -- plot as log spaced frequencies
            f1 = [freqs,[freqs(length(freqs):-1:1)]];
            f2 = [pcaproj2(1,:),[pcaproj2(2,end:-1:1)]];
            ph = fill(f1,f2,pcabkgd);hold on;
            set(gca,'FontSize',12)
            set(gca,'xtick',[10 20 40])
            xlim([freqs(1) freqs(end)])
            if cp < length(g.comps)
                set(gca,'xticklabel',[]);
            end
        end;
        set(ph,'edgecolor',pcabkgd);%set(ph,'edgealpha',0);
        if pl == (row-1)*col+2
            xlabel('Frequency (Hz)'); ylabel('Power (dB)');
        elseif pl > (row-1)*col+1
            xlabel('Frequency (Hz)');
        end;
        minl2 = min(pcaproj2(:));
        maxl2 = max(pcaproj2(:));
        if minl2 < minl
            minl = minl2;
        end;
        if maxl2 > maxl
            maxl = maxl2;
        end;
        clear pcaproj pcaproj2
    end; % to conditional plotting
    
    
    %%%%%%%%%%%%%%  Plot IM back-projections   %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    for tpp = 1:length(g.factors) % plot number of IMs for each IC
        tp = g.factors(tpp);
        backproj = winv(g.dataport,tp)*activations(tp,:); % backproj curr IM
        onebkprj = backproj(:,length(freqvec)*(rcp-1)+1:length(freqvec)*rcp);
        
        % remove outlier data by removing spectra of timewindows that exceed predefined percentage (i.e. 0.95 % of the data)
        [tpsort sortidx] = sort(winv(g.dataport,tp)); % sort timewindows by curr IM
        onebkprj = onebkprj(sortidx,:);
        onebkprj =  onebkprj + repmat(g.meanspec(rcp,:),[size(onebkprj,1) 1]);
        if g.plotperc < 1 % to shave off outer timewindows
            cutends = round(length(tpsort)-(length(tpsort)*g.plotperc));
            tpsort(1:cutends) = [];onebkprj(1:cutends,:) = [];
            tpsort(end-cutends:end) = [];onebkprj(end-cutends:end,:) = [];
        end;
        
        %% check if plotting upper/lower boundary of data, envelope or average
        if strcmp(g.plotenv, 'upper') % only pos proj
            envdata2(1,:) = onebkprj(end,:);
        elseif strcmp(g.plotenv,'lower')% only neg proj
            envdata2(1,:) = onebkprj(1,:);
        elseif strcmp(g.plotenv,'env') % plot envelope
            envdata2 = [onebkprj(end,:);onebkprj(1,:)];
        elseif strcmp(g.plotenv,'aver') % plot average
            envdata2 = mean(onebkprj,1);
        end;
        
        % If only one IM is plotted plot upper part of envelope red and
        % lower part blue
        if length(g.factors) < 2 & strcmp(g.plotenv,'env')% do red/blue if only one
            snglcols = {'r','b'};
            % for e = 1:size(envdata2,1)
            if strcmp(freqscale,'linear')
                ph = plot(freqs,envdata2(1,fr),'r-','linewidth',lnwdth); hold on;
                ph12 = gcf;
                ph12.Children(1).Children(1).Tag = taglist{1};
                hold on
                ph = plot(freqs,envdata2(2,fr),'b-','linewidth',lnwdth); hold on;
                ph12 = gcf;
                ph12.Children(1).Children(1).Tag = taglist{2};
            elseif strcmp(freqscale,'log')
                ph = semilogx(freqs,envdata2(1,fr)', 'LineWidth', 2,'Color','r');hold on
                ph12 = gcf;
                ph12.Children(1).Children(1).Tag = taglist{1};
                hold on
                ph = semilogx(freqs,envdata2(2,fr)', 'LineWidth', 2,'Color','b');hold on
                ph12 = gcf;
                ph12.Children(1).Children(1).Tag = taglist{2};
                set(gca,'FontSize',12)
                set(gca,'xtick',[10 20 40 80 120])
                xlim([freqs(1) freqs(end)])
            end;
            %                 set(ph,'color',snglcols{e});
            %                  ph12 = gcf;
            %                  ph12.Children(1).Children(1).Tag = taglist{1};
            %                  ph12.Children(2).Children(2).Tag = taglist{2};
            %                 %legendInfo{tpp} = ['IM ' num2str(g.factors(tpp))];
            %end;
        else   %% plot all IMs
            if strcmp(freqscale,'linear') % plot in linear frequency scale
                ph1 = plot(freqs,envdata2(:,fr),'k-','linewidth',2); hold on;
            elseif strcmp(freqscale,'log') % plot in log frequency scale
                ph1 = semilogx(freqs,envdata2(:,fr)', 'LineWidth', 2,'Color','m');hold on
                set(gca,'FontSize',12)
                set(gca,'xtick',[10 20 40 80 120])
                xlim([freqs(1) freqs(end)])
            end;
            set(ph1,'color',cols(tpp,:));
            ph12 = gcf;
            ph12.Children(1).Children(1).Tag = taglist{tpp};
            legendInfo{tpp} = ['IM ' num2str(g.factors(tpp))];
        end;
        collmaxes(1,tpp) = max(max(envdata2(:,fr)));
        collmins(1,tpp) = min(min(envdata2(:,fr)));
        clear plotprj envdata2
    end;
    
    %% plot mean IC spectral power
    if strcmp(freqscale,'linear') % plot in linear freq scale
        ph1 = plot(freqs,g.meanspec(rcp,fr),'k-','linewidth',lnwdth+.1); hold on;
        ph12 = gcf;
        % ph12.Children(1).Children(1).Tag = taglist{3};
    elseif strcmp(freqscale,'log') % plot in log freq scale
        ph1 = semilogx(freqs,g.meanspec(rcp,fr)', 'LineWidth', 2,'Color','k');hold on
        ph12 = gcf;
        % ph12.Children(1).Children(1).Tag = taglist{3};
        set(gca,'FontSize',12)
        set(gca,'xtick',[10 20 40 80 120])
        xlim([freqs(1) freqs(end)])
    end;
    if length(g.factors) < 2 & strcmp(g.plotenv,'env')
        ph12 = gcf;
        ph12.Children(1).Children(1).Tag = taglist{3}; % prepare figure legend
    else
        ph12 = gcf;
        ph12.Children(1).Children(1).Tag = taglist{tpp+1}; % prepare figure legend
    end
    
    if length(g.factors)> 1 && length(g.factors)< 7;
        legendInfo{tpp+1} = ['mean IC spectrum'];
        if pl <=3;
            hb = [];
            for loki = 1:tpp+1
                hb = [hb findobj(gcf,'tag',taglist{loki})];
            end
            hlegend = legend(hb,legendInfo); % plot figure legend
            set(hlegend, 'Location', 'best');
        end
    end
    
    if length(g.factors)== 1
        %ph12.Children(1).Children(1).Tag = taglist{tpp+1}; % prepare figure legend
        if strcmp(g.plotenv, 'upper') % only pos proj
            legendInfo{1} = ['IM ' num2str(g.factors(tpp)) ' upper'];
            legendInfo{2} = ['mean IC spectrum'];
        elseif strcmp(g.plotenv,'lower')
            legendInfo{1} = ['IM ' num2str(g.factors(tpp)) ' lower'];
            legendInfo{2} = ['mean IC spectrum'];
        elseif strcmp(g.plotenv,'env') % plot envelope
            legendInfo{1} = ['IM ' num2str(g.factors(tpp)) ' upper'];
            legendInfo{2} = ['IM ' num2str(g.factors(tpp)) ' lower'];
            legendInfo{3} = ['mean IC spectrum'];
        end
        if pl <=3;
            hb = [];
            for loki = 1:3
                hb = [hb findobj(gcf,'tag',taglist{loki})];
            end
            hlegend = legend(hb,legendInfo); % plot figure legend
            set(hlegend, 'Location', 'best');
        end
    end
    
    % set min max values of plot
    minl2 = min(collmins); % for use without full back-proj
    maxl2 = max(collmaxes);
    if minl2 < minl
        minl = minl2;
    end;
    if maxl2 > maxl
        maxl = maxl2;
    end;
    maxl = maxl + maxl*.02;
    minl = minl - minl*.02;
    set(gca,'ylim',[minl maxl]); set(gca,'box','off');
    
    minl = 10*(floor(minl/10)); % get to the closest 10
    maxl = 10*(ceil(maxl/10)); % get to the closest 10
    
    set(gca,'ytick',[minl:10:maxl]);        set(gca,'yticklabel',[minl:10:maxl]);
    set(gca,'ticklength',[.02 .02]);
    set(gca,'fontsize', 12)
    
    fprintf('\nComp %s of %s done.',int2str(cp),int2str(length(g.comps)));
    pl = pl+2;
    if ~isempty(g.setminmax) % if min max plotting values are provided cut plot to these values
        ylim([g.setminmax])
    end
    
end;
if col >4
    set(gcf,'Position',[100 300 1400 900]);
end
icadefs;
set(gcf,'PaperOrientation','landscape');  set(gcf,'PaperPosition',[0.25 0.25 10.5 8]);
set(gcf,'color',BACKCOLOR);
axcopy

