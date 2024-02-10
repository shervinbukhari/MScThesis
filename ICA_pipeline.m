%% Sleemory Closed-loop TMR Project: loc_mktrial
% created by H.-V.V. Ngo
%
% preprocess individual localiser runs for mvpa
%
% 1. load and filter data
% 2. segment data
% 3. visual inspection and rejection of bad trials and channels
% 4. ica for blink artifacts (on 1-Hz high pass filtered data!)
% 5. interpolate bad channels
% 6. export data


clear
ft_defaults


%% timekeeping
statime = tic;


%% directories
switch getenv('Computername')
    case 'DESKTOP-KEKPSTU'
        addpath('E:\GitHub\matlabfuns');
        
        dirRoot     = 'E:\work_uob\sleemory_tmr\data_sleep';
        dirSlSt     = 'E:\work_uob\sleemoryTMR\data_ staging';
        filMont     = 'E:\GoogleDrive\work_uob\Bham-ChanMontage\Bham-64CH-Topo.mat';
end


%% load layout
inMont = load(filMont);


%% general parameters
numSub  = 3;
timepad = [-1 3];


%% subject list and selection
inSub       = cell(numSub,2);
inSub(1,:)  = {'sleemory04'; 's04_hg'};
inSub(2,:)  = {'sleemory05'; 's05_ps'};
inSub(3,:)  = {'sleemory09'; 's09_ap'};


iSub = input(sprintf('Choose subject (n = %d) to preprocess: ',numSub));

fprintf('Process subject %d/%d: %s\n', iSub, numSub, inSub{iSub,1});


%% load data
% ... unfiltered data
fprintf(' load unfiltered data\n');
tfg                     = [];
tfg.dataset             = fullfile(dirData,[inSub{iSub} '_sleep.vhdr']);
tfg.channel             = {'all'; '-E*'; '-M*'};
goddata_raw             = ft_preprocessing(tfg);

%... filter data
fprintf(' create filtered dataset\n');
tfg                     = [];
tfg.hpfilter            = 'yes';
tfg.hpfreq              = 1;
tfg.hpfiltord           = 3;
tfg.hpinstabilityfix    = 'reduce';
tfg.lpfilter            = 'yes';
tfg.lpfreq              = 35;
tfg.lpfiltord           = 5;
tfg.lpinstabilityfix    = 'reduce';
tfg.bsfilter            = 'yes';
tfg.bsfreq              = [48 52];
goddata_filt            = ft_preprocessing(tfg,goddata_raw);

%... useful parameters
fsample = goddata_raw.fsample;
datalen = goddata_raw.sampleinfo(2);


%% import sleep staging
tfg             = [];
tfg.slStLen     = 30;
tfg.hypnogram   = fullfile(dirRoot,[inSub{iSub} '_staging.txt']);
tfg.type        = 'SchlafAus';

staging = hvn_importHypnogram(tfg, fsample, datalen);

    
%% segmentation
stagefltr = ismember(staging,[5 15]);        %% get REM sleep stage only
sleepbouts = hvn_extrctBnryBouts(stagefltr);                          %% Determine sleepbouts

%... prepare first coarse-grain segmentation
tfg             = [];
tfg.minlength   = 8.196;
tfg.trl         = [sleepbouts(:,1), sleepbouts(:,2), zeros(size(sleepbouts,1),1)];
        
trl_raw     = ft_redefinetrial(tfg,goddata_raw);
trl_filt    = ft_redefinetrial(tfg,goddata_filt);


%... further segmentation in fine-grained segments
tfg         = [];
tfg.length  = 8.196;
tfg.overlap = 0.5;

trl_raw     = ft_redefinetrial(tfg,trl_raw);
trl_filt    = ft_redefinetrial(tfg,trl_filt);

%... add trial number for later matching between raw and filtered data
trl_raw.trialinfo   = cat(2,trl_raw.trialinfo,(1:size(trl_raw.trialinfo,1))');
trl_filt.trialinfo  = trl_raw.trialinfo;
    

%% artifact inspection and rejection on filtered data
%... save trial number and channel labels
trl_old = (1 : numMrkr)';
ch_old  = trl_filt.label;

%... summary
tfg         = [];
tfg.layout  = inMont.lay;
tfg.method  = 'summary';
tfg.alim    = 250;
trl_filt    = ft_rejectvisual(tfg, trl_filt);
    
%... trials
tfg         = [];
tfg.layout  = inMont.lay;
tfg.method  = 'trial';
tfg.alim    = 250;
trl_filt    = ft_rejectvisual(tfg, trl_filt);

%... save bad trials and bad channels
bad_ch  = setdiff(ch_old, trl_filt.label);
bad_trl = setdiff(trl_old,trl_filt.trialinfo(:,3));


%% remove bad trials and channels in unfiltered data
tfg         = [];
tfg.channel = trl_filt.label;
tfg.trials  = trl_filt.trialinfo(:,3);

trl_raw = ft_preprocessing(tfg,trl_raw);


%% ICA for eog artifact removal
%... calculate ICA
tfg         = [];
tfg.method  = 'runica';
all_comp    = ft_componentanalysis(tfg, trl_filt);

%... save components
save(fullfile(dirRoot,[inSub{iSub} '_ICA.mat']),'-struct','all_comp');

%... plot components as topographies
figure;
tfg             = [];
tfg.zlim        = 'maxabs';
tfg.component   = 1:30;         % specify the component(s) that should be plotted
tfg.layout      = inMont.lay;
tfg.comment     = 'no';

ft_topoplotIC(tfg,all_comp)

figure;
tfg             = [];
tfg.zlim        = 'maxabs';
tfg.component   = 31 : numel(all_comp.label);
tfg.layout      = inMont.lay;
tfg.comment     = 'no';

ft_topoplotIC(tfg, all_comp)


%... plot component times series
tfg             = [];
tfg.viewmode    = 'component';
tfg.layout      = inMont.lay;
tfg.continuous  = 'no';
tfg.blocksize   = 4;
tfg.channel     = 1:numel(all_comp.label);
ft_databrowser(tfg, all_comp);

clear trl_filt

pause


%... user input to remove components
rmvComp = input('Specfiy components to be removed (leave empty for none!): ');

%... removal of unwanted ica components if necessary
if ~isempty(rmvComp)
    tfg             = [];
    tfg.component   = rmvComp; % vector with excl. components

    trl_raw = ft_rejectcomponent(tfg, all_comp, trl_raw);
end

clear allcomp
   

%% interpolate bad channels if necessary
if ~isempty(bad_ch)
    tfg                 = [];
    tfg.method          = 'spline';
    tfg.missingchannel  = bad_ch;
    tfg.neighbours      = inMont.neighbours;
    tfg.layout          = inMont.lay;

    trl_raw = ft_channelrepair(tfg, trl_raw);
end

    
%% save results
%... eeg data
save(fullfile(dirRoot,[inSub{iSub} '_trls.mat']),'-struct', 'trl_raw');

%... supplementary
outSupplmt          = [];
outSupplmt.rmvComp  = rmvComp;
outSupplmt.bad_ch   = bad_ch;
outSupplmt.bad_trl  = bad_trl;

save(fullfile(dirRoot,[inSub{iSub} '_supplmt.mat']),'-struct', 'outSupplmt');

clear outSupplmt trl_raw


%% print summary
fprintf('and to summarize (n = %d):\n', numel(bad_ch));
fprintf('removed the following channels: ');
for iCh = 1 : numel(bad_ch)
    fprintf('%s ', bad_ch{iCh});
end
fprintf('\n');
fprintf('and the trials (n = %d): ', numel(bad_trl));
for iTrl = 1 : numel(bad_trl)
    fprintf('%d ', bad_trl(iTrl));
end
fprintf('\n');


%% Timekeeping
fprintf('shit went down in %.2f s\n', toc(statime));
