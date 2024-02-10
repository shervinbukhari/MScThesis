
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

%% REDEFINE WITH OVERLAP
for N=1:numel(Data_fltr_altica)
cfg = [];
cfg.length =5;
cfg.overlap =0.5;
Ovlp{N}=ft_redefinetrial(cfg, Data_fltr_altica{N});
end

cfg =[];
cfg.keepsampleinfo  = 'yes';
SVLDs2 = ft_appenddata(cfg, Data_fltr_altica{1,5}, Data_fltr_altica{1,6}, Data_fltr_altica{1,7});
notLDs2 = ft_appenddata(cfg, Data_fltr_altica{1,1}, Data_fltr_altica{1,2}, Data_fltr_altica{1,3}, Data_fltr_altica{1,4});

%% Power Spectrum Analysis
for N=1:numel(Ovlp)
cfg.method = 'mtmfft';
cfg.keeptrials = 'yes';                            
cfg.output      = 'powandcsd';
cfg.taper = 'hanning';
cfg.channel ={'all' '-H*' '-V*'};
cfg.foi = 0:(1/5):45;
cfg.pad     = 'nextpow2';
cfg.tapsmofrq = 0;
orig{N}=ft_freqanalysis(cfg, Ovlp{N});
cfg.method='irasa';
frac=ft_freqanalysis(cfg, Ovlp{N});
end
%only SVLDs
cfg.method = 'mtmfft';
cfg.keeptrials = 'yes';                            
cfg.output      = 'powandcsd';
cfg.taper = 'dpss';
cfg.channel ={'all' '-H*' '-V*'};
cfg.foi = 0:(1/5):45;
cfg.tapsmofrq = 2;
tfg.pad     = 'nextpow2';
SVLD_orig=ft_freqanalysis(cfg, SVLDs2);
%cfg.method = 'irasa';
%SVLD_frac=ft_freqanalysis(cfg, SVLDs2);
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
notLD_orig=ft_freqanalysis(cfg, notLDs2);
% cfg.method = 'irasa';
% notLD_frac=ft_freqanalysis(cfg, notLDs2);
%% subtract the fractal component from the power spectrum
for N=1:numel(orig)
cfg               = [];
cfg.parameter     = 'powspctrm';
cfg.operation     = 'x2-x1';
osci{N} = ft_math(cfg,frac{N}, orig{N});
end

cfg               = [];
cfg.parameter     = 'powspctrm';
cfg.operation     = 'x2-x1';
SVLD_osci = ft_math(cfg, SVLD_frac, SLVD_orig);

cfg               = [];
cfg.parameter     = 'powspctrm';
cfg.operation     = 'x2-x1';
notLD_osci = ft_math(cfg, notLD_frac, notLD_orig);

%% Oscillatory Power of Original for each frequency band
%SVLD
cfg = [];
cfg.frequency = [1 4];
orig_SVLD_d=ft_selectdata(cfg, SVLD_orig);
cfg = [];
cfg.frequency = [4 8];
orig_SVLD_t=ft_selectdata(cfg, SVLD_orig);
cfg = [];
cfg.frequency =[8 12];
orig_SVLD_a=ft_selectdata(cfg, SVLD_orig);
cfg = [];
cfg.frequency =[12 20];
orig_SVLD_b=ft_selectdata(cfg, SVLD_orig);
cfg = [];
cfg.frequency =[20 45];
orig_SVLD_g=ft_selectdata(cfg, SVLD_orig);
%nonLD
cfg = [];
cfg.frequency = [1 4];
orig_notLD_d=ft_selectdata(cfg, notLD_orig);
cfg = [];
cfg.frequency = [4 8];
orig_notLD_t=ft_selectdata(cfg, notLD_orig);
cfg = [];
cfg.frequency =[8 12];
orig_notLD_a=ft_selectdata(cfg, notLD_orig);
cfg = [];
cfg.frequency =[12 20];
orig_notLD_b=ft_selectdata(cfg, notLD_orig);
cfg = [];
cfg.frequency =[20 45];
orig_notLD_g=ft_selectdata(cfg, notLD_orig);
%% Oscillatory Power for each frequency band
%SVLD
cfg = [];
cfg.frequency = [1 4];
osci_SVLD_d=ft_selectdata(cfg, SVLD_osci);
cfg = [];
cfg.frequency = [4 8];
osci_SVLD_t=ft_selectdata(cfg, SVLD_osci);
cfg = [];
cfg.frequency =[8 12];
osci_SVLD_a=ft_selectdata(cfg, SVLD_osci);
cfg = [];
cfg.frequency =[12 20];
PSD_SVLD_b=ft_selectdata(cfg, SVLD_osci);
cfg = [];
cfg.frequency =[20 45];
PSD_SVLD_g=ft_selectdata(cfg, SVLD_osci);
%nonLD
cfg = [];
cfg.frequency = [1 4];
PSD_notLD_d=ft_selectdata(cfg, notLD_osci);
cfg = [];
cfg.frequency = [4 8];
PSD_notLD_t=ft_selectdata(cfg, notLD_osci);
cfg = [];
cfg.frequency =[8 12];
PSD_notLD_a=ft_selectdata(cfg, notLD_osci);
cfg = [];
cfg.frequency =[12 20];
PSD_notLD_b=ft_selectdata(cfg, notLD_osci);
cfg = [];
cfg.frequency =[20 45];
PSD_notLD_g=ft_selectdata(cfg, notLD_osci);

%%
for N=1:numel(orig)
figure;
cfg = [];
% cfg.showlabels   = 'yes';
cfg.layout = lay;
cfg.colorbar     = 'yes';
xlabel('Frequency (Hz)');
ylabel('Total power (µV^2)');
ft_singleplotER(cfg, orig{N})
title('PSD')
end
%%
%Descriptives
for N=1:numel(orig)
cfg =[];
PSD_dsc{N}=ft_freqdescriptives(cfg, orig{N});
end 
cfg =[];
cfg.keeptrials    = 'yes';
cfg.variance      = 'yes';
PSD_notLD_dsc = ft_freqdescriptives(cfg, not);
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
ylabel('Absolute power (µV^2)');
title('Average PSD of all SVLD trials')
ft_topoplotER(cfg, PSD_SVLD_a)

% finally, single subject lines colored by group and error bars (SEM) and stat differences
fh = figure;
data1 = SVLD_orig;
data2 = notLD_orig;
f= linspace(0,45,226);

mp1 = mean(SVLD_orig.powspctrm,[1 2]);
var1 = var(mean(SVLD_orig.powspctrm,2));
var1 = squeeze(var1);
n1 = 112;
sem1 = sqrt(var1)./sqrt(n1);
mp1=squeeze(mp1)
mp2 = mean(notLD_orig.powspctrm,[1 2]);
mp2=squeeze(mp2)
var2 = var(mean(notLD_orig.powspctrm,2));
var2=squeeze(var2)
n1 = 112;
sem2 = sqrt(var2)./sqrt(n1);
figure;
% ft_plot_line(f,mp1,'color',[0.5 0.5 1],'linewidth',1.5)
% 
% ft_plot_line(f,mp2,'color',[1 0.5 0.5],'linewidth',1.5)


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
mp3=mean(LD_orig_contrast.powspctrm, [1 2]);
mp3=squeeze(mp3);
var3=var(mean(LD_orig_contrast.powspctrm,2));
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
cfg.ylim =[-0.25 0.55];
ft_singleplotER(cfg, LD_orig_contrast)
%% Make contrast structures for topoplotting
%Average PSD SVLD > non-lucid
LD_orig_contrast = SVLD_orig;
LD_orig_contrast.powspctrm = SVLD_orig.powspctrm - notLD_orig.powspctrm;
LD_orig_contrast_a = orig_SVLD_a;
LD_orig_contrast_a.powspctrm = orig_SVLD_a.powspctrm - orig_notLD_a.powspctrm;
LD_orig_contrast_d = orig_SVLD_d;
LD_orig_contrast_d.powspctrm = orig_SVLD_d.powspctrm - orig_notLD_d.powspctrm;
LD_orig_contrast_b = orig_SVLD_b;
LD_orig_contrast_b.powspctrm = orig_SVLD_b.powspctrm - orig_notLD_b.powspctrm;
LD_orig_contrast_t = orig_SVLD_t;
LD_orig_contrast_t.powspctrm = orig_SVLD_t.powspctrm - orig_notLD_t.powspctrm;
LD_orig_contrast_g = orig_SVLD_g;
LD_orig_contrast_g.powspctrm = orig_SVLD_g.powspctrm - orig_notLD_g.powspctrm;


%%
cfg = [];
cfg.layout = lay;
cfg.colorbar     = 'SouthOutside';
% cfg.zlim = [-0.0025 0.0025]
cfg.marker = 'off'
cfg.markersymbol ='.';
cfg.markersize = 8;
cfg.comment          = 'no';
% figure('position',[680 240 1039 420]);
subplot(1,5,1); ft_topoplotER(cfg, LD_orig_contrast_d); title('Delta (1-4 Hz)');
ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
colormap(flipud(brewermap(64,'RdBu'))) % change the colormap
subplot(1,5,2); ft_topoplotER(cfg, LD_orig_contrast_t); title('Theta (4-8 Hz)');
ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
colormap(flipud(brewermap(64,'RdBu'))) % change the colormap
subplot(1,5,3); ft_topoplotER(cfg, LD_orig_contrast_a); title('Alpha (8-12 Hz)');
ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
colormap(flipud(brewermap(64,'RdBu'))) % change the colormap
cfg.marker ='off';
subplot(1,5,4); ft_topoplotER(cfg, LD_orig_contrast_b); title('Beta (12-20 Hz)');
ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
colormap(flipud(brewermap(64,'RdBu'))) % change the colormap
subplot(1,5,5); ft_topoplotER(cfg, LD_orig_contrast_g); title('Gamma (20-45 Hz)');
ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
colormap(flipud(brewermap(64,'RdBu'))) % change the colormap

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
stats_delta_orig = ft_freqstatistics(cfg, orig_SVLD_d, orig_notLD_d);
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
% cfg.highlightsizeseries = ['6' '6']
cfg.layout = lay;
cfg.colobar='yes';
ft_clusterplot(cfg, stats_delta_orig);title('Delta (1-4 Hz)'); hold on;
ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
colormap(flipud(brewermap(64,'RdBu'))) % change the colormap
ft_clusterplot(cfg,stats_theta_orig);title ('Theta (4-8 Hz)');
ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
colormap(flipud(brewermap(64,'RdBu'))) % change the colormap
ft_clusterplot(cfg,stats_alpha_orig);title('Alpha (8-12 Hz)');
ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
colormap(flipud(brewermap(64,'RdBu'))) % change the colormap
ft_clusterplot(cfg, stats_beta_orig);title('Beta (12-20 Hz)');
ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
colormap(flipud(brewermap(64,'RdBu'))) % change the colormap
ft_clusterplot(cfg, stats_gamma_orig);title('Gamma (20-45 Hz)');
ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
colormap(flipud(brewermap(64,'RdBu'))) % change the colormap
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

