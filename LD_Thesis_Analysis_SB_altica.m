
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
%load('LD_data') %Load Data
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
for N=1:numel(Data_fltr_altica)
cfg = [];
cfg.length =5;
cfg.overlap =0.5;
Data_ovlp2{N}=ft_redefinetrial(cfg, Data_fltr_altica{N});
end

cfg =[];
cfg.keepsampleinfo  = 'yes'
SVLDs2 = ft_appenddata(cfg, Data_fltr_altica{1,5}, Data_fltr_altica{1,6}, Data_fltr_altica{1,7});
notLDs2 = ft_appenddata(cfg, Data_fltr_altica{1,1}, Data_fltr_altica{1,2}, Data_fltr_altica{1,3}, Data_fltr_altica{1,4})

%% Take first derivative of the data for 1/f normalisation
% load('LD_Data')
%1/f 5s
for N=1:numel(Data_ovlp2)
cfg=[]
cfg.derivative = 'yes'
Data_deriv2{N}=ft_preprocessing(cfg, Data_ovlp2{N});
end

cfg = [];
cfg.length =5;
cfg.overlap =0.5;
SVLD_ovlp2=ft_redefinetrial(cfg, SVLDs2);

cfg =[];
cfg.length =5;
cfg.overlap =0.5;
notLDs_ovlp2=ft_redefinetrial(cfg, notLDs2);

%1/f 1s
for N=1:numel(Data_fltr)
cfg =[]
cfg.derivative = 'yes'
tfa_clean_deriv{N}=ft_preprocessing(cfg, Data_fltr{N})
end

% SVLD
cfg=[]
cfg.derivative = 'yes'
SVLD_deriv2=ft_preprocessing(cfg, SVLD_ovlp2)

%notLD
cfg=[]
cfg.derivative = 'yes'
notLD_deriv2=ft_preprocessing(cfg, notLDs_ovlp2)

%% Power Spectrum Analysis
for N=1:numel(Data_deriv2)
cfg.method = 'mtmfft';
cfg.keeptrials = 'yes';                            
cfg.output      = 'powandcsd';
cfg.taper = 'dpss';
cfg.channel ={'all' '-H*' '-V*'};
cfg.foi = 0:(1/5):45;
cfg.pad     = 'nextpow2';
cfg.tapsmofrq = 2;
Data_PSD2{N}=ft_freqanalysis(cfg, Data_deriv2{N});
end

%only SVLDs
cfg.method = 'mtmfft';
cfg.keeptrials = 'yes';                            
cfg.output      = 'powandcsd';
cfg.taper = 'dpss';
cfg.channel ={'all' '-H*' '-V*'};
cfg.foi = 0:(1/5):45;
cfg.tapsmofrq = 2;
cfg.pad     = 'nextpow2';
PSD_SVLD2=ft_freqanalysis(cfg, SVLD_deriv2);

%notLD
cfg.method = 'mtmfft';
cfg.keeptrials = 'yes';                            
cfg.output      = 'powandcsd';
cfg.trials = randperm(120,108);
cfg.taper = 'dpss';
cfg.channel ={'all' '-H*' '-V*'};
cfg.foi = 0:(1/5):45;
cfg.pad     = 'nextpow2';
cfg.tapsmofrq = 2;
PSD_notLD2=ft_freqanalysis(cfg, notLD_deriv2);

%% For each frequency band
%SVLD
cfg = [];
cfg.frequency = [1 4];
PSD_SVLD_d=ft_selectdata(cfg, PSD_SVLD2);
cfg = [];
cfg.frequency = [4 8];
PSD_SVLD_t=ft_selectdata(cfg, PSD_SVLD2);
cfg = [];
cfg.frequency =[8 12];
PSD_SVLD_a=ft_selectdata(cfg, PSD_SVLD2);
cfg = [];
cfg.frequency =[12 20];
PSD_SVLD_b=ft_selectdata(cfg, PSD_SVLD2);
cfg = [];
cfg.frequency =[20 45];
PSD_SVLD_g=ft_selectdata(cfg, PSD_SVLD2);
%nonLD
cfg = [];
cfg.frequency = [1 4];
PSD_notLD_d=ft_selectdata(cfg, PSD_notLD2);
cfg = [];
cfg.frequency = [4 8];
PSD_notLD_t=ft_selectdata(cfg, PSD_notLD2);
cfg = [];
cfg.frequency =[8 12];
PSD_notLD_a=ft_selectdata(cfg, PSD_notLD2);
cfg = [];
cfg.frequency =[12 20];
PSD_notLD_b=ft_selectdata(cfg, PSD_notLD2);
cfg = [];
cfg.frequency =[20 45];
PSD_notLD_g=ft_selectdata(cfg, PSD_notLD2);

%%
%PSD for each nap avg over time and channels
for N=1:numel(Data_PSD2)
figure;
cfg = [];
% cfg.showlabels   = 'yes';
cfg.layout = lay;
cfg.colorbar     = 'yes';
xlabel('Frequency (Hz)');
ylabel('normalized power (µV^2)');
ft_topoplotER(cfg, Data_PSD2{N})
title('PSD')
end
%%
%Descriptives
for N=1:numel(Data_PSD2)
cfg =[];
PSD_dsc{N}=ft_freqdescriptives(cfg, Data_PSD2{N});
end 
cfg =[];
cfg.keeptrials    = 'yes';
cfg.variance      = 'yes';
PSD_notLD_dsc = ft_freqdescriptives(cfg, PSD_notLD2);
PSD_SVLD_dsc = ft_freqdescriptives(cfg, PSD_SVLD2);
%% Grand average
cfg=[];
PSD_dsc_noLD2=ft_freqgrandaverage(cfg, PSD_dsc2{1,1},PSD_dsc2{1,2},PSD_dsc2{1,3},PSD_dsc2{1,4});
PSD_dsc_SVLD2=ft_freqgrandaverage(cfg, PSD_dsc2{1,5}, PSD_dsc2{1,6},PSD_dsc2{1,7})
%% PSD SINGLEPLOT ER 
cfg = []
cfg.showlabels   = 'yes';
cfg.layout = lay;
cfg.colorbar     = 'yes';
cfg.linewidth = 1.5
figure;
xlabel('Frequency (Hz)');
ylabel('normalized power (µV^2)');
title('Average Power Spectral Density')
ft_singleplotER(cfg, PSD_notLD2, PSD_SVLD2);
legend('non-LDs','SVLDs')
ft_singleplotER(cfg, PSD_LD_contrast2)
hold on;
%Plot SVLD
cfg = [];
cfg.showlabels   = 'yes';
cfg.layout = lay;
figure;
xlabel('Frequency (Hz)');
ylabel('Normalized power (µV^2)');
title('Average PSD of all SVLD trials')
ft_topoplotER(cfg, PSD_SVLD_a)

% finally, single subject lines colored by group and error bars (SEM) and stat differences
fh = figure;
data1 = PSD_SVLD2;
data2 = PSD_notLD2;
f= linspace(0,45,226)

mp1 = mean(PSD_SVLD2.powspctrm,[1 2]);
var1 = var(mean(PSD_SVLD2.powspctrm,2));
var1 = squeeze(var1)
n1 = 112;
sem1 = sqrt(var1)./sqrt(n1);
mp1=squeeze(mp1)
mp2 = mean(PSD_notLD2.powspctrm,[1 2]);
mp2=squeeze(mp2)
var2 = var(mean(PSD_notLD2.powspctrm,2));
var2=squeeze(var2)
n1 = 112;
sem2 = sqrt(var2)./sqrt(n1);
figure;
ft_plot_line(f,mp1,'color',[0.5 0.5 1],'linewidth',1.5)

ft_plot_line(f,mp2,'color',[1 0.5 0.5],'linewidth',1.5)


x = [f flip(f)];
facealpha = 0.4;

y_above = mp1+sem1;
y_below = mp1-sem1;
y_patch = [y_above; flip(y_below)];
ft_plot_patch(x, y_patch,'facecolor','r','facealpha',facealpha)
y_above = mp2+sem2;
y_below = mp2-sem2;
y_patch = [y_above; flip(y_below)];
ft_plot_patch(x, y_patch,'facecolor','b','facealpha',facealpha)

p1=ft_plot_line(f, mp1,'color','r','linewidth',2);
p2=ft_plot_line(f, mp2,'color','b','linewidth',2);
legend([p1 p2],{'SVLDs','non-LDs'})
title('Average Power Spectral Density')
mp3=mean(PSD_LD_contrast2.powspctrm, [1 2]);
mp3=squeeze(mp3);
var3=var(mean(PSD_LD_contrast2.powspctrm,2));
var3=squeeze(var3);
sem3=sqrt(var3)./sqrt(n1);
facealpha = 0.50;

y_above = mp3+sem3;
y_below = mp3-sem3;
y_patch = [y_above; flip(y_below)];
ft_plot_patch(x, y_patch,'facecolor','k','facealpha',facealpha)
cfg = [];
cfg.linewidth = 1.5;
cfg.linecolor = 'k';
xlabel('Frequency (Hz)');
ylabel('Normalized power (µV^2)');
title('Power Spectral Density Difference')
cfg = [];
cfg.linewidth = 1.5;
cfg.linecolor = 'k';
cfg.ylim =[-0.0004 0.0004];
ft_singleplotER(cfg, PSD_LD_contrast2)
%% Make contrast structures for topoplotting
%Average PSD SVLD > non-lucid
PSD_LD_contrast2 = PSD_SVLD2;
PSD_LD_contrast2.powspctrm = PSD_SVLD2.powspctrm - PSD_notLD2.powspctrm;
PSD_LD_contrast_a = PSD_SVLD_a;
PSD_LD_contrast_a.powspctrm = PSD_SVLD_a.powspctrm - PSD_notLD_a.powspctrm;
PSD_LD_contrast_d = PSD_SVLD_d;
PSD_LD_contrast_d.powspctrm = PSD_SVLD_d.powspctrm - PSD_notLD_d.powspctrm;
PSD_LD_contrast_b = PSD_SVLD_b;
PSD_LD_contrast_b.powspctrm = PSD_SVLD_b.powspctrm - PSD_notLD_b.powspctrm;
PSD_LD_contrast_t = PSD_SVLD_t;
PSD_LD_contrast_t.powspctrm = PSD_SVLD_t.powspctrm - PSD_notLD_t.powspctrm;
PSD_LD_contrast_g = PSD_SVLD_g;
PSD_LD_contrast_g.powspctrm = PSD_SVLD_g.powspctrm - PSD_notLD_g.powspctrm;


%%
cfg = [];
cfg.layout = lay;
cfg.colorbar     = 'SouthOutside';
cfg.zlim = [-0.0025 0.0025]
cfg.marker = 'on'
cfg.markersymbol ='.';
cfg.markersize = 8;
cfg.comment          = 'no';
% figure('position',[680 240 1039 420]);
subplot(1,5,1); ft_topoplotER(cfg, PSD_LD_contrast_d); title('Delta (1-4 Hz)');
ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
colormap(flipud(brewermap(64,'RdBu'))) % change the colormap
subplot(1,5,2); ft_topoplotER(cfg, PSD_LD_contrast_t); title('Theta (4-8 Hz)');
ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
colormap(flipud(brewermap(64,'RdBu'))) % change the colormap
subplot(1,5,3); ft_topoplotER(cfg, PSD_LD_contrast_a); title('Alpha (8-12 Hz)');
ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
colormap(flipud(brewermap(64,'RdBu'))) % change the colormap
cfg.marker ='off';
subplot(1,5,4); ft_topoplotER(cfg, PSD_LD_contrast_b); title('Beta (12-20 Hz)');
ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
colormap(flipud(brewermap(64,'RdBu'))) % change the colormap
subplot(1,5,5); ft_topoplotER(cfg, PSD_LD_contrast_g); title('Gamma (20-45 Hz)');
ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
colormap(flipud(brewermap(64,'RdBu'))) % change the colormap
xlabel('Frequency (Hz)');
ylabel('normalized power (µV^2)');
title('Power Difference')
%%
for N=1:numel(Data_PSD2)
figure;
cfg = [];
cfg.interactive = 'yes';
cfg.layout = lay;
xlabel('Frequency (Hz)');
ylabel('normalized power (µV^2)');
ft_singleplotER(cfg, Data_PSD2{N})
end

%% Permutation Testing
% do NOT EXECUTE this yet, it is just to introduce the function´
cfg                 = [];
cfg.channel          = 'all';
cfg.avgoverchan     = 'no';
cfg.frequency        = 'all';
cfg.avgoverfreq     = 'yes';
cfg.parameter        = 'powspctrm';
cfg.method           = 'ft_statistics_montecarlo';
cfg.statistic        = 'ft_statfun_depsamplesT'; %SHOULD BE DEPENDENT T 
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.minnbchan        = 2;
cfg.tail             = 0;
cfg.clustertail      = cfg.tail;
cfg.neighbours = neighbours;
cfg.alpha            = 0.05;
cfg.correcttail      = 'alpha';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.computeprob      = 'yes';
cfg.numrandomization = 1000;
design = zeros(1,size(PSD_SVLD2.powspctrm,1) + size(PSD_notLD2.powspctrm,1))
design(1,1:size(PSD_SVLD2.powspctrm,1)) = 1; % 1=LUCID %0=NLD
design(1,(size(PSD_SVLD2.powspctrm,1)+1):(size(PSD_SVLD2.powspctrm,1)+...
size(PSD_notLD2.powspctrm,1))) = 2;
 design(2,:)=1:216;
design(2,1:end)=1;
cfg.design=design;
cfg.ivar = 1;
cfg.uvar = 2;
% cfg.wvar = within session trials
stats_gamma = ft_freqstatistics(cfg, PSD_SVLD_g, PSD_notLD_g);
%% Cluster plot
cfg = [];
PSD_SVLD_freqdcrpt = ft_freqdescriptives(cfg, PSD_SVLD);
PSD_notLD_freqdcrpt  = ft_freqdescriptives(cfg, PSD_notLD);
stats.raweffect = PSD_SVLD_freqdcrpt.powspctrm - PSD_notLD_freqdcrpt.powspctrm;
PSD_LD_contrast1 = PSD_SVLD;
PSD_LD_contrast1.powspctrm = PSD_SVLD.powspctrm - PSD_noLD.powspctrm;
stats_alpha.raweffect = PSD_LD_contrast_a.powspctrm;

cfg = [];
cfg.alpha  = 0.05;
cfg.interactive = 'yes';
cfg.highlightcolorneg = [1 1 1];
cfg.highlightsymbolseries = ['.' '+'];
cfg.layout = lay;
cfg.colobar='yes';
ft_clusterplot(cfg, stats_delta);title('Delta (1-4 Hz)'); hold on;
ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
colormap(flipud(brewermap(64,'RdBu'))) % change the colormap
ft_clusterplot(cfg,stats_theta);title ('Theta (4-8 Hz)');
ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
colormap(flipud(brewermap(64,'RdBu'))) % change the colormap
ft_clusterplot(cfg,stats_alpha);title('Alpha (8-12 Hz)');
ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
colormap(flipud(brewermap(64,'RdBu'))) % change the colormap
ft_clusterplot(cfg, stats_beta);title('Beta (12-20 Hz)');
ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
colormap(flipud(brewermap(64,'RdBu'))) % change the colormap
ft_clusterplot(cfg, stats_gamma);title('Gamma (20-45 Hz)');
ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
colormap(flipud(brewermap(64,'RdBu'))) % change the colormap
%% Cluster Analyis Hongi
stats_alpha.powspctrm = PSD_LD_contrast_a.powspctrm;

cfg = [];
cfg.parameter ='stat';
cfg.highlight ='on';
cfg.highlightsymbol ='*';
cfg.highlightschannels = stats_alpha.mask;
cfg.layout = lay;
ft_topoplotER(cfg, stats_alpha)

%% Source Analysis
load('standard_bem'),load('elec_aligned_bem'),load('leadfield_bem')
elec1005=ft_read_sens('standard_1005.elc');
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

