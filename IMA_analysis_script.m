% Test cases Independent Modulators analysis (IMA)

% addpath('/Volumes/Drive/IMAT_project/eeglab2019_1/')
% addpath('/Volumes/Drive/IMAT_project/IMAT_plugin/mfiles/IMAT/')

% Launch EEGLAB
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

% Load sample data
filepath = '/Users/amon-ra/WORK/Data/IMAT_data/';
filename = 'RestEC_S03_ContAMICAdip.set';

EEG = pop_loadset('filename',filename,'filepath',filepath);
[ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
eeglab redraw;

%% Single subject analysis
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Running IMA
[EEG, IMA] = pop_runIMA(EEG, 'freqscale', 'linear',...
                             'frqlim', [6 120],...
                             'pcfac', 8,...
                             'cycles', [6 0.5],...
                             'selectICs', {'brain'},...
                             'icatype', 'infomax');

%% Plotting Options here
%--------------------------------------------------------------------------

%% First plot and options: pop_plotspecdecomp
pop_plotspecdecomp(EEG, 'plottype', 'comb', 'comps', [1 2 3 4 5 6]) % need to check the numbers of ICs to plot when empty
pop_plotspecdecomp(EEG, 'plottype', 'ics', 'comps', [1 2 3 4 5 6])
pop_plotspecdecomp(EEG, 'plottype', 'ims', 'comps', [1 2 3 4 5 6])

%% Second plot: pop_plotspecenv
pop_plotspecenv(EEG,'comps', [1 2 3 4 5 6], 'factors', [1 2 4 5], 'frqlim', [3 40], 'plotenv', 'upper');%, 'setminmax', [30 80]);

%% Third plot: pop_plotIMtimecourse
pop_plotIMtimecourse(EEG, 'comps', [1 2 3], 'factors', [1 5], 'frqlim', [6 40], 'smoothing', 0.1) % maybe put plots on one panel

%'plotcond', 'off',
%   collect_templates_clustering(IMA)
%  collect_templates_clustering_Templatesel(IMA)

%% Multiple subjects analysis(STUDY)
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%% Loading STUDY
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
[STUDY ALLEEG] = pop_loadstudy('filename', 'S3_RestECEO.study', 'filepath', '/Users/amon-ra/WORK/Data/IMAT_data');
CURRENTSTUDY = 1; EEG = ALLEEG; CURRENTSET = [1:length(EEG)];
eeglab redraw;
%%
% Running IMA
[STUDY] = pop_runIMA_study(STUDY, ALLEEG, 'freqscale', 'log',...
                                          'frqlim', [6 120],...
                                          'pcfac', 8,...
                                          'cycles', [6 0.5],...
                                          'selectICs', {'brain'},...
                                          'icatype', 'infomax');

% Plotting Options here
%--------------------------------------------------------------------------
%% Plot 1: pop_plotspecdecomp_study
pop_plotspecdecomp_study(STUDY, 'plottype', 'comb')

%% Plot 3: pop_plotspecenv_study
pop_plotspecenv_study(STUDY,'comps', [1 3 5 6], 'factors', [1 2 3 5 7], 'frqlim', [3 40],'plotcond', 'on');

%% Plot 4: pop_plotIMtimecourse_study
pop_plotIMtimecourse_study(STUDY, 'comps', [1 3 5], 'factors', [1 5], 'frqlim', [6 40], 'plotcond', 'on', 'smoothing', 0.1,...
    'plotIMtf', 'on', 'plotIMtime', 'on')


%   peakrange = [8 13];
%   stretch_spectra = 'on';
%   targetpeakfreq = 10;

%% Clustering 
% [TempKeep, TempKeepInd, dipsources, scalpmaps] = collect_templates_clustering_Templatesel(IMA, peakrange, stretch_spectra, targetpeakfreq);
[IMA] = pop_collecttemplates(STUDY, 'peakrange', [8 14],...
                                    'stretch_spectra', 'on',...
                                    'targetpeakfreq', 10,...
                                    'plot_templ', 'on');
%'peakrange', [9 13], '
[STUDY] = pop_clusterIMAtemplates(STUDY, ALLEEG, 'nclust', 3);

[STUDY] = pop_subclusterIMAtemplates(STUDY, 'clust', [2 3], 'nclust', 2);

%% Plot 2: pop_plotIMAcluster
pop_plotIMAcluster(STUDY, 'plotclust','on', 'plottemplates', 'on', 'plotscalpmaps', 'on', 'plotdipsources', 'on', 'plotsubclusters', 'on')




