
%% INITIALIZATION
clear all; close all;

%%% PLEASE MAKE SURE:
%%% 1. YOU HAVE Matlab2013b or later
%%% 2. YOU HAVE ALL REQUIRED TOOLBOXES (typically the case)
%%% 3. YOUR SLEEPTRIP FOLDER DIRECTLY HAS THE FILES IN THEM, 
%%%    e.g. while extracting a zip from github maybe you need to refer to 
%%%    D:/sleeptrip-master/sleeptrip-master instead of D:/sleeptrip-master
%%% 4. NO FieldTrip NOR related toolboxes like EEGLAB are loaded or in the 
%%%    PATH of your MATLAB

pathToFieldTrip = 'M:\Documents\FieldTrip\fieldtrip';

addpath(pathToFieldTrip);
addpath('P:\3013077.01\');
addpath('P:\3013077.01\Scripts_301307701\Shervin')
load('LD_data') %Load Data
load('P:\3013077.01\Scripts_301307701\Shervin\topography\cue_topo.mat') %Load Topo
%%% disable some toolboxes if necessary and check if 
%%% 'signal_toolbox', 'signal_blocks' are available 
%%% because they are helpful to have.
%  toggleToolbox('names');
%  toggleToolbox('dsp','off');
%  toggleToolbox('all','query');
%  license('inuse');
ft_defaults
%%
cfg = [];
cfg.viewmode = 'vertical';
cfg.ylim = [-50 50]
cfg.channel = [1:10]
cfg.channelclamped = ['HEOG1','VEOG1'];
ft_databrowser(cfg, Data_fltr_eyes{1,1});

%% REDEFINE WITH OVERLAP
for N=1:numel(Data_fltr)
cfg = [];
cfg.length =5;
cfg.overlap =0.5;
Data_ovlp{N}=ft_redefinetrial(cfg, Data_fltr{N});
end

cfg =[]
cfg.keepsampleinfo  = 'yes'
SVLDs = ft_appenddata(cfg, Data_fltr{1,5}, Data_fltr{1,6}, Data_fltr{1,7});
notLDs = ft_appenddata(cfg, Data_fltr{1,1}, Data_fltr{1,2}, Data_fltr{1,3}, Data_fltr{1,4})

%% Take first derivative of the data for 1/f normalisation
% load('LD_Data')
%1/f 5s
for N=1:numel(Data_ovlp)
cfg=[]
cfg.derivative = 'yes'
Data_deriv{N}=ft_preprocessing(cfg, Data_ovlp{N});
end

cfg = [];
cfg.length =5;
cfg.overlap =0.5;
SVLD_ovlp=ft_redefinetrial(cfg, SVLDs);

cfg =[];
cfg.length =5;
cfg.overlap =0.5;
notLDs_ovlp=ft_redefinetrial(cfg, notLDs);

%1/f 1s
for N=1:numel(Data_fltr)
cfg =[]
cfg.derivative = 'yes'
tfa_clean_deriv{N}=ft_preprocessing(cfg, Data_fltr{N})
end

% SVLD
cfg=[]
cfg.derivative = 'yes'
SVLD_deriv=ft_preprocessing(cfg, SVLD_ovlp)

%notLD
cfg=[]
cfg.derivative = 'yes'
notLD_deriv=ft_preprocessing(cfg, notLDs_ovlp)

%% Power Spectrum Analysis
for N=1:numel(Data_deriv)
cfg.method = 'mtmfft';
cfg.keeptrials = 'yes';                            
cfg.output      = 'powandcsd';
cfg.taper = 'dpss'
cfg.channel ={'all' '-H*' '-V*'};
cfg.foi = 1:1:45;
cfg.tapsmofrq = 0.5;
cfg.output = 'powandcsd'
Data_PSD{N}=ft_freqanalysis(cfg, Data_deriv{N});
end

%only SVLDs
cfg.method = 'mtmfft';
cfg.keeptrials = 'yes';                            
cfg.output      = 'pow';
cfg.taper = 'dpss';
cfg.channel ={'all' '-H*' '-V*'};
cfg.foi = 1:1:45
cfg.tapsmofrq = 0.5;
cfg.output = 'powandcsd'
PSD_SVLD=ft_freqanalysis(cfg, SVLD_deriv);

%notLD
cfg.method = 'mtmfft';
cfg.keeptrials = 'yes';                            
cfg.output      = 'pow';
cfg.taper = 'dpss'
cfg.channel ={'all' '-H*' '-V*'};
cfg.foi = 1:1:45
cfg.tapsmofrq = 0.5;
cfg.output = 'powandcsd'
PSD_notLD=ft_freqanalysis(cfg, notLD_deriv);

%% For each frequency band
%Delta
cfg.method = 'mtmfft';
cfg.keeptrials = 'yes';                            
cfg.output      = 'pow';
cfg.taper = 'dpss';
cfg.channel ={'all' '-H*' '-V*'};
cfg.foi = 1:1:4
cfg.tapsmofrq = 0.5;
PSD_SVLD_d=ft_freqanalysis(cfg, SVLD_deriv);

cfg.method = 'mtmfft';
cfg.keeptrials = 'yes';                            
cfg.output      = 'pow';
cfg.taper = 'dpss';
cfg.channel ={'all' '-H*' '-V*'};
cfg.foi = 1:1:4
cfg.tapsmofrq = 0.5;
PSD_notLD_d=ft_freqanalysis(cfg, notLD_deriv);


%Theta
cfg.method = 'mtmfft';
cfg.keeptrials = 'yes';                            
cfg.output      = 'pow';
cfg.taper = 'dpss';
cfg.channel ={'all' '-H*' '-V*'};
cfg.foi = 4:1:8
cfg.tapsmofrq = 0.5;
PSD_SVLD_t=ft_freqanalysis(cfg, SVLD_deriv);

cfg.method = 'mtmfft';
cfg.keeptrials = 'yes';                            
cfg.output      = 'pow';
cfg.taper = 'dpss';
cfg.channel ={'all' '-H*' '-V*'};
cfg.foi = 4:1:8
cfg.tapsmofrq = 0.5;
PSD_notLD_t=ft_freqanalysis(cfg, notLD_deriv);

%Alpha
cfg.method = 'mtmfft';
cfg.keeptrials = 'yes';                            
cfg.output      = 'pow';
cfg.taper = 'dpss';
cfg.channel ={'all' '-H*' '-V*'};
cfg.foi = 8:1:12
cfg.tapsmofrq = 0.5;
PSD_notLD_a=ft_freqanalysis(cfg, notLD_deriv);

cfg.method = 'mtmfft';
cfg.keeptrials = 'yes';                            
cfg.output      = 'pow';
cfg.taper = 'dpss';
cfg.channel ={'all' '-H*' '-V*'};
cfg.foi = 8:1:12
cfg.tapsmofrq = 0.5;
PSD_SVLD_a=ft_freqanalysis(cfg, SVLD_deriv);

%Beta1
cfg.method = 'mtmfft';
cfg.keeptrials = 'yes';                            
cfg.output      = 'pow';
cfg.taper = 'dpss';
cfg.channel ={'all' '-H*' '-V*'};
cfg.foi = 12:1:20
cfg.tapsmofrq = 0.5;
PSD_SVLD_b=ft_freqanalysis(cfg, SVLD_deriv);

cfg.method = 'mtmfft';
cfg.keeptrials = 'yes';                            
cfg.output      = 'pow';
cfg.taper = 'dpss';
cfg.channel ={'all' '-H*' '-V*'};
cfg.foi = 12:1:20
cfg.tapsmofrq = 0.5;
PSD_notLD_b=ft_freqanalysis(cfg, notLD_deriv);

%Gamma
cfg.method = 'mtmfft';
cfg.keeptrials = 'yes';                            
cfg.output      = 'pow';
cfg.taper = 'dpss';
cfg.channel ={'all' '-H*' '-V*'};
cfg.foi = 20:1:45
cfg.tapsmofrq = 0.5;
PSD_SVLD_g=ft_freqanalysis(cfg, SVLD_deriv);

cfg.method = 'mtmfft';
cfg.keeptrials = 'yes';                            
cfg.output      = 'pow';
cfg.taper = 'dpss';
cfg.channel ={'all' '-H*' '-V*'};
cfg.foi = 20:1:45
cfg.tapsmofrq = 0.5;
PSD_notLD_g=ft_freqanalysis(cfg, notLD_deriv);

%%
% save(...)
load('P:\3013077.01\Scripts_301307701\Shervin\Data Files\Freq_data.mat')

%PSD for each nap avg over time and channels
for N=1:numel(Data_PSD)
figure;
cfg = [];
% cfg.showlabels   = 'yes';
cfg.layout = lay;
cfg.colorbar     = 'yes';
xlabel('Frequency (Hz)');
ylabel('absolute power (µV^2)');
ft_singleplotER(cfg, Data_PSD{N})
title('PSD')
end
%% Raw & Normalized Absolute Power
%Average PSD over each exp
for N=1:numel(Data_PSD)
cfg =[];
PSD_avg{N}=ft_freqdescriptives(cfg, Data_PSD{N})
end 

cfg=[];
PSD_avg_noLD=ft_freqgrandaverage(cfg, PSD_avg{1,1},PSD_avg{1,2},PSD_avg{1,3},PSD_avg{1,4});
PSD_avg_SVLD=ft_freqgrandaverage(cfg, PSD_avg{1,5}, PSD_avg{1,6},PSD_avg{1,7})
PSD_avg_TWC =ft_freqgrandaverage(cfg, PSD_avg{1,5}, PSD_avg{1,6}, PSD_avg{1,8})
PSD_avg_LDWTC =ft_freqgrandaverage(cfg, PSD_avg{1,5}, PSD_avg{1,6});

cfg = [];
cfg.showlabels   = 'yes';
cfg.layout = lay;
cfg.colorbar     = 'yes';
cfg.linewidth = 1.5
cfg.graphcolor = 'k';
figure;
xlabel('Frequency (Hz)');
ylabel('absolute power (µV^2)');
title('Average Power Spectral Density')
ft_singleplotER(cfg, PSD_avg_noLD, PSD_avg_SVLD);
legend('non-LDs','SVLDs')
ft_singleplotER(cfg, PSD_LD_contrast1)

%Plot SVLD
cfg = [];
cfg.showlabels   = 'yes';
cfg.layout = lay;
figure;
xlabel('Frequency (Hz)');
ylabel('absolute power (µV^2)');
title('Average PSD of all SVLD trials')
ft_topoplotER(cfg, PSD_avg_SVLD)

%Average PSD SVLD > non-lucid
PSD_LD_contrast1 = PSD_avg_SVLD;
PSD_LD_contrast1.powspctrm = PSD_avg_SVLD.powspctrm - PSD_avg_noLD.powspctrm;

cfg = [];
cfg.showlabels   = 'yes';
cfg.interactive = 'yes'
cfg.lay = lay;
cfg.colorbar     = 'yes';
cfg.highlightsymbol  = {'o','*'};
cfg.highlightcolor   = [0 0 0];
cfg.highlightsize    = 6;
cfg.markersymbol     = '.';
cfg.comment          = 'no';
figure;
xlabel('Frequency (Hz)');
ylabel('absolute power (µV^2)');
title('SVLDs > non-LDs')
ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
colormap(flipud(brewermap(64,'RdBu'))) % change the colormap
ft_movieplotER(cfg, PSD_LD_contrast1)

for N=1:numel(Data_PSD)
figure;
cfg = [];
cfg.showlabels   = 'yes';
cfg.interactive = 'yes'
cfg.layout = lay;
cfg.colorbar     = 'yes';
cfg.highlightsymbol  = {'o','*'};
cfg.highlightcolor   = [0 0 0];
cfg.highlightsize    = 6;
cfg.markersymbol     = '.';
cfg.comment          = 'no';
cfg.colormap         = 'jet';
cfg.colorbar         = 'yes'
figure;
xlabel('Frequency (Hz)');
ylabel('absolute power (µV^2)');
title('PSD SVLD > non-LDs')
ft_topoplotER(cfg, Data_PSD{N})
end

%% Permutation Testing
% do NOT EXECUTE this yet, it is just to introduce the function´
cfg                 = [];
cfg.channel          = 'all';
cfg.avgoverchan     = 'no';
cfg.frequency        = 'all';
cfg.avgoverfreq     = 'no';
cfg.parameter        = 'powspctrm';
cfg.method           = 'ft_statistics_montecarlo';
cfg.statistic        = 'ft_statfun_indepsamplesT'; %SHOULD BE DEPENDENT T 
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.clusterthreshold = 'nonparametric_common';
cfg.minnbchan        = 2;
cfg.tail             = 0;
cfg.clustertail      = cfg.tail;
cfg.neighbours = neighbours;
cfg.alpha            = 0.05;
cfg.correcttail      = 'alpha';
cfg.computeprob      = 'yes';
cfg.numrandomization = 100000;
design = zeros(1,size(PSD_SVLD.powspctrm,1) + size(PSD_notLD.powspctrm,1))
design(1,1:size(PSD_SVLD.powspctrm,1)) = 1; % 1=LUCID %0=NLD
design(1,(size(PSD_SVLD.powspctrm,1)+1):(size(PSD_SVLD.powspctrm,1)+...
size(PSD_notLD.powspctrm,1))) = 2;
% design(2,:)=1:108;
% design(2,109:end)=1
cfg.design=design;
cfg.ivar = 1;
% cfg.uvar = 2;
% cfg.wvar = within session trials
stats = ft_freqstatistics(cfg, PSD_SVLD, PSD_notLD);
%save freq_stats stats

%% Cluster plot
cfg = [];
PSD_SVLD_freqdcrpt = ft_freqdescriptives(cfg, PSD_SVLD);
PSD_notLD_freqdcrpt  = ft_freqdescriptives(cfg, PSD_notLD);
stats.raweffect = PSD_SVLD_freqdcrpt.powspctrm - PSD_notLD_freqdcrpt.powspctrm;
PSD_LD_contrast1 = PSD_avg_SVLD;
PSD_LD_contrast1.powspctrm = PSD_avg_SVLD.powspctrm - PSD_avg_noLD.powspctrm;
stats_alpha.raweffect = PSD_avg_SVLD.powspctrm - PSD_avg_noLD.powspctrm;

cfg = [];
cfg.alpha  = 0.05;
cfg.parameter = 'raweffect';
cfg.interactive = 'yes';
cfg.layout = lay;
ft_clusterplot(cfg, stats_alpha);
% 
%% Source Analysis
load('standard_bem'),load('elec_aligned_bem'),load('leadfield_bem')
elec1005=ft_read_sens('standard_1005.elc')
ft_plot_mesh(vol.bnd(1,:)); hold on; ft_plot_sens(elec_aligned_bem3, 'label','label')

ft_plot_mesh(vol.bnd(1,:)); hold on; ft_plot_sens(elec,'orientation','true')

cfg = [];
cfg.headmodel= vol;
cfg.elec = elec1005;
cfg.reducerank = 3;
leadfield_bem = ft_prepare_leadfield(cfg);

% create spatial filter using the lcmv beamformer
cfg                  = [];
cfg.method           = 'dics';
cfg.sourcemodel      = leadfield_bem; % leadfield
cfg.elec = elec1005;
cfg.headmodel        = vol; % volume conduction model (headmodel)
cfg.lcmv.keepfilter  = 'yes';
cfg.lcmv.fixedori    = 'yes'; % project on axis of most variance using SVD
source               = ft_sourceanalysis(cfg, PSD_SVLD);
%%
%% PERMUTE CROSS SPECTRUM
cfg                 = [];
cfg.channel          = 'all';
cfg.avgovergchan     = 'no';
cfg.frequency        = 'all';
cfg.avgoverfreq     = 'no';
cfg.parameter        = 'crsspctrm';
cfg.method           = 'ft_statistics_montecarlo';
cfg.statistic        = 'ft_statfun_indepsamplesT'; %SHOULD BE DEPENDENT T 
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.clusterthreshold = 'nonparametric_common';
cfg.minnbchan        = 2;
cfg.tail             = 0;
cfg.clustertail      = cfg.tail;
cfg.neighbours = neighbours;
cfg.alpha            = 0.05;
cfg.correcttail      = 'alpha';
cfg.computeprob      = 'yes';
cfg.numrandomization = 5000;
design = zeros(1,size(PSD_SVLD.powspctrm,1) + size(PSD_notLD.powspctrm,1))
design(1,1:size(PSD_SVLD.powspctrm,1)) = 1; % 1=LUCID %0=NLD
design(1,(size(PSD_SVLD.powspctrm,1)+1):(size(PSD_SVLD.powspctrm,1)+...
size(PSD_notLD.powspctrm,1))) = 2;
% design(2,:)=1;
cfg.design=design;
cfg.ivar = 1;
% cfg.uvar = 2;
stats_crosspsctrm = ft_freqstatistics(cfg, PSD_SVLD, PSD_notLD);

