
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

%% Power Spectrum Analysis
%only SVLDs
cfg.method = 'mtmfft';
cfg.keeptrials = 'yes';                            
cfg.output      = 'pow';
cfg.taper = 'hanning';
cfg.channel ={'all' '-H*' '-V*'};
cfg.foi = 0:(1/5):30;
cfg.tapsmofrq = 0;
cfg.pad     = 'nextpow2';
SVLD_hann=ft_freqanalysis(cfg, SVLD_ovlp2);
cfg.foi=30:(1/5):45;
cfg.taper ='dpss';
cfg.tapsmofrq=4;
SVLD_dpss=ft_freqanalysis(cfg, SVLD_ovlp2);

%notLD
cfg.method = 'mtmfft';
cfg.keeptrials = 'yes';                            
cfg.output      = 'pow';
cfg.trials = randperm(120,108);
cfg.taper = 'hann';
cfg.channel ={'all' '-H*' '-V*'};
cfg.foi = 0:(1/5):30;
cfg.pad     = 'nextpow2';
cfg.tapsmofrq = 0;
notLD_hann=ft_freqanalysis(cfg, notLDs_ovlp2);
cfg.foi=30:(1/5):45;
cfg.taper ='dpss';
cfg.tapsmofrq=4;
notLD_dpss=ft_freqanalysis(cfg, notLDs_ovlp2);
%smooth Hann
notLD_hann.powspctrm = smoothdata(notLD_hann.powspctrm,3,'movmean',3);
SVLD_hann.powspctrm = smoothdata(SVLD_hann.powspctrm,3,'movmean',3);
%append
cfg=[];
cfg.parameter ='powspctrm';
cfg.appenddim ='freq';
SVLD_hs=ft_appendfreq(cfg, SVLD_hann, SVLD_dpss);
notLD_hs=ft_appendfreq(cfg, notLD_hann, notLD_dpss);

%% For each frequency band
%SVLD
cfg = [];
cfg.frequency = [1 4];
SVLD_hs_d=ft_selectdata(cfg, SVLD_hs);
cfg = [];
cfg.frequency = [4 8];
SVLD_hs_t=ft_selectdata(cfg, SVLD_hs);
cfg = [];
cfg.frequency =[8 12];
SVLD_hs_a=ft_selectdata(cfg, SVLD_hs);
cfg = [];
cfg.frequency =[12 30];
SVLD_hs_b=ft_selectdata(cfg, SVLD_hs);
cfg = [];
cfg.frequency =[30 45];
SVLD_hs_g=ft_selectdata(cfg, SVLD_hs);
%nonLD
cfg = [];
cfg.frequency = [1 4];
notLD_hs_d=ft_selectdata(cfg, notLD_hs);
cfg = [];
cfg.frequency = [4 8];
notLD_hs_t=ft_selectdata(cfg, notLD_hs);
cfg = [];
cfg.frequency =[8 12];
notLD_hs_a=ft_selectdata(cfg, notLD_hs);
cfg = [];
cfg.frequency =[12 30];
notLD_hs_b=ft_selectdata(cfg, notLD_hs);
cfg = [];
cfg.frequency =[30 45];
notLD_hs_g=ft_selectdata(cfg, notLD_hs);
%% PSD SINGLEPLOT ER 
% finally, single subject lines colored by group and error bars (SEM) and stat differences

f= linspace(0,45,226);
SVLD_hs.logpowspctrm =log10(SVLD_hs.powspctrm);
notLD_hs.logpowspctrm=log10(notLD_hs.powspctrm);

mp1 = mean(SVLD_hs.logpowspctrm,[1 2]);
var1 = var(mean(SVLD_hs.logpowspctrm,2));
var1 = squeeze(var1);
n1 = 112;
sem1 = sqrt(var1)./sqrt(n1);
mp1=squeeze(mp1);
mp2 = mean(notLD_hs.logpowspctrm,[1 2]);
mp2=squeeze(mp2);
var2 = var(mean(notLD_hs.logpowspctrm,2));
var2=squeeze(var2);
n1 = 112;
sem2 = sqrt(var2)./sqrt(n1);
figure;

x = [f flip(f)];
facealpha = 0.4;

y_above = mp1+sem1;
y_below = mp1-sem1;
y_patch = [y_above; flip(y_below)];
figure;ft_plot_patch(x, y_patch,'facecolor','r','facealpha',facealpha)
y_above = mp2+sem2;
y_below = mp2-sem2;
y_patch = [y_above; flip(y_below)];
ft_plot_patch(x, y_patch,'facecolor','b','facealpha',facealpha)

p1=ft_plot_line(f, mp1,'color','r','linewidth',2);
p2=ft_plot_line(f, mp2,'color','b','linewidth',2);
legend([p1 p2],{'SVLDs','nonLDs'})
title('Average Power Spectral Density')

mp3=mean(LD_contrast_hs.logpowspctrm, [1 2]);
mp3=squeeze(mp3);
var3=var(mean(LD_contrast_hs.logpowspctrm,2));
var3=squeeze(var3);
sem3=sqrt(var3)./sqrt(n1);
facealpha = 0.50;
y_above = mp3+sem3;
y_below = mp3-sem3;
y_patch = [y_above; flip(y_below)];
figure;ft_plot_patch(x, y_patch,'facecolor','k','facealpha',facealpha)
ft_plot_line(f, mp3,'color','k','linewidth',2);
cfg = [];
cfg.linewidth = 1.5;
cfg.linecolor = 'k';
xlabel('Frequency (Hz)');
ylabel('Normalized power (µV^2)');
title('Power Spectral Density Difference')
cfg = [];
cfg.linewidth = 1.5;
cfg.linecolor = 'k';
cfg.ylim =[-0.5 2];
ft_singleplotER(cfg, PSD_LD_contrast_hs)
ylabel('Log Power (µV^2)')
title('Power Spectral Density')
%% Make contrast structures for topoplotting
%Average PSD SVLD > non-lucid
LD_contrast_hs = SVLD_hs;
LD_contrast_hs.powspctrm = SVLD_hs.powspctrm - notLD_hs.powspctrm;
LD_contrast_hs.logpowspctrm=SVLD_hs.logpowspctrm - notLD_hs.logpowspctrm;
LD_contrast_a_hs = SVLD_hs_a;
LD_contrast_a_hs.powspctrm = SVLD_hs_a.powspctrm - notLD_hs_a.powspctrm;
LD_contrast_d_hs = SVLD_hs_d;
LD_contrast_d_hs.powspctrm = SVLD_hs_d.powspctrm - notLD_hs_d.powspctrm;
LD_contrast_b_hs = SVLD_hs_b;
LD_contrast_b_hs.powspctrm = SVLD_hs_b.powspctrm - notLD_hs_b.powspctrm;
LD_contrast_t_hs = SVLD_hs_t;
LD_contrast_t_hs.powspctrm = SVLD_hs_t.powspctrm - notLD_hs_t.powspctrm;
LD_contrast_g_hs = SVLD_hs_g;
LD_contrast_g_hs.powspctrm = SVLD_hs_g.powspctrm - notLD_hs_g.powspctrm;

%%
cfg = [];
cfg.layout = lay;
cfg.zlim = 'maxabs';
cfg.colorbar     = 'SouthOutside';
cfg.marker = 'off';
cfg.markersymbol ='.';
cfg.markersize = 8;
cfg.comment          = 'no';
figure;subplot(1,5,1); ft_topoplotER(cfg, LD_contrast_d_hs); title('Delta (1-4 Hz)');
ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
colormap(flipud(brewermap(64,'RdBu'))) % change the colormap
subplot(1,5,2); ft_topoplotER(cfg, LD_contrast_t_hs); title('Theta (4-8 Hz)');
ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
colormap(flipud(brewermap(64,'RdBu'))) % change the colormap
subplot(1,5,3); ft_topoplotER(cfg, LD_contrast_a_hs); title('Alpha (8-12 Hz)');
ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
colormap(flipud(brewermap(64,'RdBu'))) % change the colormap
cfg.marker ='off';
subplot(1,5,4); ft_topoplotER(cfg, LD_contrast_b_hs); title('Beta (12-30 Hz)');
ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
colormap(flipud(brewermap(64,'RdBu'))) % change the colormap
subplot(1,5,5); ft_topoplotER(cfg, LD_contrast_g_hs); title('Gamma (30-45 Hz)');
ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
colormap(flipud(brewermap(64,'RdBu'))) % change the colormap
%%
cfg =[];
figure;ft_singleplotER(cfg, LD_contrast_d_hs);
figure;ft_singleplotER(cfg, SVLD_hs_d, notLD_hs_d);
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
cfg.minnbchan        = 4;
cfg.tail             = 0;
cfg.clustertail      = cfg.tail;
cfg.neighbours = neighbours;
cfg.alpha            = 0.05;
cfg.correcttail      = 'alpha';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'wcm';
cfg.computeprob      = 'yes';
cfg.numrandomization = 1000;
design = zeros(1,size(SVLD_hs.powspctrm,1) + size(notLD_hs.powspctrm,1));
design(1,1:size(SVLD_hs.powspctrm,1)) = 1; % 1=LUCID %0=NLD
design(1,(size(SVLD_hs.powspctrm,1)+1):(size(SVLD_hs.powspctrm,1)+...
size(notLD_hs.powspctrm,1))) = 2;
 design(2,:)=1:216;
design(2,1:end)=1;
cfg.design=design;
cfg.ivar = 1;
cfg.uvar = 2;
% cfg.wvar = within session trials
st_hs_g_min3 = ft_freqstatistics(cfg, SVLD_hs_g, notLD_hs_g);
%% Cluster plot
cfg = [];
cfg.alpha=0.05;
cfg.layout = lay;
cfg.highlightseries = {'labels' 'off' 'off' 'off' 'off'};
ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
colormap(flipud(brewermap(64,'RdBu'))) % change the colormap
figure;ft_clusterplot(cfg, st_hs_a);
%% Per band
cfg = [];
cfg.alpha  = 0.05;
cfg.interactive = 'yes';
cfg.zlim = 'maxabs';
cfg.highlightcolorneg = [1 1 1];
cfg.highlightsymbolseries = ['.' 'o'];
cfg.highlightsizeseries =[6 3];
cfg.layout = lay;
cfg.colobar='SouthOutside';
ft_clusterplot(cfg, st_hs_d);title('Delta (1-4 Hz)'); hold on;
ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
colormap(flipud(brewermap(64,'RdBu'))) % change the colormap
ft_clusterplot(cfg,st_hs_t);title ('Theta (4-8 Hz)');
ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
colormap(flipud(brewermap(64,'RdBu'))) % change the colormap
ft_clusterplot(cfg,st_hs_a);title('Alpha (8-12 Hz)');
ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
colormap(flipud(brewermap(64,'RdBu'))) % change the colormap
ft_clusterplot(cfg, st_hs_b);title('Beta (12-20 Hz)');
ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
colormap(flipud(brewermap(64,'RdBu'))) % change the colormap
ft_clusterplot(cfg, st_hs_g);title('Gamma (20-45 Hz)');
ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
colormap(flipud(brewermap(64,'RdBu'))) % change the colormap
