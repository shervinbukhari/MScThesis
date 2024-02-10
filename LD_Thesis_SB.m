%%                      Master Script for LD Analysis
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


%%% add sleeptrip to the path
pathToSleepTrip = 'M:\Documents\Sleeptrip\sleeptrip-master\sleeptrip-master';
pathToFieldTrip = 'M:\Documents\FieldTrip\fieldtrip';

addpath(pathToSleepTrip);
addpath(pathToFieldTrip);
addpath('P:\3013077.01\');
addpath('P:\3013077.01\Scripts_301307701\Shervin')

%%% disable some toolboxes if necessary and check if 
%%% 'signal_toolbox', 'signal_blocks' are available 
%%% because they are helpful to have.
%  toggleToolbox('names');
%  toggleToolbox('dsp','off');
%  toggleToolbox('all','query');
%  license('inuse');
ft_defaults

%% LOAD DATA, ADD INFO AND DEFINE REM PERIODS AND OFFSET
sr=500; %sampling rate
sub = [];
sub.name = '25LDES';
sub.data = [];
sub.data.dataset = {
    ['P:\3013077.01\pp6_25LDES\25LDES_EEG_nap1\25LDES_nap1_Change Sampling Rate.eeg']...
    ['P:\3013077.01\pp6_25LDES\25LDES_EEG_nap1\25LDES_nap1.2_Change Sampling Rate.eeg'] ['P:\3013077.01\pp6_25LDES\25LDES_EEG_nap1\25LDES_nap1.2_Change Sampling Rate.eeg'] ['P:\3013077.01\pp6_25LDES\25LDES_EEG_nap1\25LDES_nap1.2_Change Sampling Rate.eeg']...
    ['P:\3013077.01\pp6_25LDES\25LDES_EEG_DC2\25LDES_nap2_DC2.eeg'] ['P:\3013077.01\pp6_25LDES\25LDES_EEG_DC2\25LDES_nap2_DC2.eeg']...
    ['P:\3013077.01\pp6_25LDES\25LDES_EEG_DC3\25LDES_nap3_DC3.eeg']...
    ['P:\3013077.01\pp6_25LDES\25LDES_DC4\25LDES_nap4.eeg'] ['P:\3013077.01\pp6_25LDES\25LDES_DC4\25LDES_nap4.eeg']};

sub.REM_on =    {(sr*870)... %REM onset in samples
              (sr*2271) (sr*3450) (sr*4404)...
              (sr*405) (sr*6921)...
              (sr*6351)...
              (sr*840) (sr*5310)};


sub.REM_off =    {(sr*1158)... %REM end in samples
               (sr*2577) (sr*3657) (sr*4553)...
               (sr*648) (sr*7872)...
               (sr*6534)...
               (sr*1049) (sr*5568)};
sub.REM_offset =    {(sr*0)... %Offset = 0
               (sr*0) (sr*0) (sr*0)...
               (sr*0) (sr*0)...
               (sr*0)...
               (sr*0),(sr*0)};
%% SPECIFY BAD CHANNELS AND MISSING CHANNELS (NEXT SECTION REQUIRED TO RUN BEFORE THIS STEP)
sub.data.badchannel{1,1} = {'C3' 'CP1' 'CP3' 'FFC6h' 'CPP4h' 'TTP8h' 'T8' 'O1' 'P10'}% 'SI3' 'SI5' 'SI6' 'SI4','IIz' 'VEOG1' 'VEOG2' 'HEOG1' 'HEOG2' 'EMG1' 'EMG2' 'EMG3' 'ECG' }; %, %'SI6', 'SI4'
sub.data.badchannel{1,2} ={'C3' 'CCP3h' 'Pz' 'T8' 'CP1' 'TTP8h' 'CPP4h' 'P10' 'O1' }%'SI3' 'SI5' 'SI6' 'SI4','IIz'  'VEOG1' 'VEOG2' 'HEOG1' 'HEOG2' 'EMG1' 'EMG2' 'EMG3' 'ECG'}'; %' 
sub.data.badchannel{1,3} = {'C3' 'O1' 'T8' 'Pz' 'CP1' 'TTP8h' 'CCP3h' 'CPP4h' 'P10' }%'SI3' 'SI5' 'SI6' 'SI4','IIz'  'VEOG1' 'VEOG2' 'HEOG1' 'HEOG2' 'EMG1' 'EMG2' 'EMG3' 'ECG'}';
sub.data.badchannel{1,4} = {'CP1' 'CP3' 'O1' 'Pz' 'CP3' 'FCC6h' 'CCP3h' 'CCP4h' 'CPP4h' 'P10' }%'SI3' 'SI5' 'SI6' 'SI4','IIz' 'VEOG1' 'VEOG2' 'HEOG1' 'HEOG2' 'EMG1' 'EMG2' 'EMG3' 'ECG'}';
sub.data.badchannel{1,5} = {}%'SI3' 'SI5' 'SI6' 'SI4','IIz'  'VEOG1' 'VEOG2' 'HEOG1' 'HEOG2' 'EMG1' 'EMG2' 'EMG3' 'ECG'};
sub.data.badchannel{1,6} ={}%'SI3' 'SI5' 'SI6' 'SI4','IIz'  'VEOG1' 'VEOG2' 'HEOG1' 'HEOG2' 'EMG1' 'EMG2' 'EMG3' 'ECG'};
sub.data.badchannel{1,7} = {}%'SI3' 'SI5' 'SI6' 'SI4','IIz'  'VEOG1' 'VEOG2' 'HEOG1' 'HEOG2' 'EMG1' 'EMG2' 'EMG3' 'ECG'};

sub.data.badchannel{1,8} ={'O1' 'T7' 'FC1' 'F5' 'C1' 'CP1' 'CP3' 'CP6' 'TP8'...
'P1' 'PO3' 'AFF5h' 'FFC1h' 'FFC2h' 'FFC2h' 'FFC3h' 'FFC5h' 'CCP3h' 'CCP4h' ...
'CCP6h' 'CPP1h' 'CPP2h' 'CPP4h' 'CPP5h' 'PPO6h' 'TPP8h' 'TTP8h' 'TTP7h'...
'PPO1h' 'PPO2h' 'P9' 'P10' 'Iz' 'SI3' 'SI5' 'SI6' 'SI4','IIz'  'VEOG1' 'VEOG2' 'HEOG1' 'HEOG2' 'EMG1' 'EMG2' 'EMG3' 'ECG'}';
sub.data.badchannel{1,9} = sub.data.badchannel{1,8}';

for NapNumber=1:numel(sub.data.badchannel)
Good_Chans{NapNumber}= setdiff(preprocs_raw{NapNumber}.label, sub.data.badchannel{NapNumber})
end;
save('Good_Chans', 'Good_Chans')
%% READING IN THE DATA AND VARIOUS PREPROCESS CONFIGURATIONS

%Just Reading, No Referencing or Processing 
for NapNumber = 1:numel(sub.data.dataset);
    cfg = [];
    cfg.dataset = sub.data.dataset{NapNumber};
    start = sub.REM_on{NapNumber};
    stop  = sub.REM_off{NapNumber};
    offset = sub.REM_offset{NapNumber};
    cfg.trl = [start stop offset];
    cfg.trl = round(cfg.trl);
    preproc = ft_preprocessing(cfg);
    preprocs_raw{NapNumber} = preproc;
end

%Seperate Ext. Elecs
for NapNumber=1:numel(preprocs_raw)
cfg = [];
cfg.trials = 'all';
cfg.channel = {'VEOG*'; 'M*'; 'HEOG*'; 'E*' 'Pz'};
ext_chans{NapNumber} = ft_selectdata(cfg,preprocs_raw{NapNumber});
end
save('25LDES_raw', 'preprocs_raw', 'ext_chans')

%Referencing EMG
for NapNumber =1:numel(ext_chans)
cfg = [];
cfg.channel = {'EMG*'};
cfg.reref = 'yes';
cfg.refchannel= 'EMG3';
cfg.hpfilter = 'yes';
cfg.hpfreq = 10;
emg_chans{NapNumber}=ft_preprocessing(cfg,ext_chans{NapNumber});
end

%Processing and Referencing HEOG
for NapNumber =1:numel(ext_chans)
cfg = [];
cfg.channel = {'H*'}
cfg.reref = 'yes';
cfg.refchannel = 'HEOG2';
cfg.hpfilter = 'yes';
cfg.lpfilter = 'yes';
cfg.hpfreq = 0.3;
cfg.lpfreq = 45;
heog_chans_nofltr{NapNumber} =ft_preprocessing(cfg, ext_chans{NapNumber}); 
end

%Processing and Referencing HEOG
for NapNumber =1:numel(ext_chans)
cfg = [];
cfg.channel = {'V*'}
cfg.reref = 'yes';
cfg.refchannel = 'VEOG2';
cfg.hpfilter = 'yes';
cfg.lpfilter = 'yes';
cfg.hpfreq = 0.3;
cfg.lpfreq = 45;
veog_chans_nofltr{NapNumber} =ft_preprocessing(cfg, ext_chans{NapNumber}); 
end

%Derive Radial EOG
for NapNumber =1:numel(ext_chans)
cfg = [];
cfg.channel = {'V*'}
cfg.reref = 'yes';
cfg.refchannel = 'VEOG2';
cfg.hpfilter = 'yes';
cfg.lpfilter = 'yes';
cfg.hpfreq = 0.3;
cfg.lpfreq = 45;
veog_chans_nofltr{NapNumber} =ft_preprocessing(cfg, ext_chans{NapNumber}); 
end
    
%Combine EOGs and Discard reference chans
for NapNumber =1:numel(eog_chans_nofltr)
    cfg = []
    cfg.channel = ['VEOG1','HEOG1'];
    eog_chans_nofltr{NapNumber} = ft_appenddata(cfg, veog_chans_nofltr{NapNumber},heog_chans_nofltr{NapNumber});
end

for NapNumber =1:numel(eog_chans_nofltr)
    cfg = []
    cfg.channel = ({'VEOG1','HEOG1'});
    eog_chans_nofltr{NapNumber} = ft_selectdata(cfg, eog_chans_nofltr{NapNumber});
end

%ECG
for NapNumber =1:numel(ext_chans)
cfg =[];
cfg.channel = 'ECG'; %no ref ch.
ecg_chans{NapNumber} = ft_preprocessing(cfg,ext_chans{NapNumber});
end
save ('25LDES_extchans_ref', 'eog_chans','emg_chans','ecg_chans')


%Processing for Channel Inspection 
for NapNumber = 1:numel(sub.data.dataset);
    cfg = [];
    cfg.dataset = sub.data.dataset{NapNumber};
    start = sub.REM_on{NapNumber};
    stop  = sub.REM_off{NapNumber};
    offset = sub.REM_offset{NapNumber};
    cfg.trl = [start stop offset];
    cfg.trl = round(cfg.trl);
    cfg.reref = 'yes';    
    cfg.refchannel = {'M1','M2'};
    cfg.bpfilter = 'yes';
    cfg.bsfilter = 'yes';
    cfg.bpfreq = [0.3 35];
    cfg.bsfreq = [48 52; 98 102];
    preproc = ft_preprocessing(cfg);
    preprocs_Chinspct{NapNumber} = preproc;
end
save('25LDES_Chinspct','preprocs_Chinspct')

%Processing for Artifact Rejection and ICA
for NapNumber = 1:numel(sub.data.dataset);
    cfg = [];
    cfg.dataset = sub.data.dataset{NapNumber};
    start = sub.REM_on{NapNumber};
    stop  = sub.REM_off{NapNumber};
    offset = sub.REM_offset{NapNumber};
    cfg.hpfilter = 'yes';
    cfg.hpfreq              = 1;
    cfg.hpfiltord           = 3;
    cfg.hpinstabilityfix    = 'reduce';
    cfg.channel = Good_Chans{NapNumber};
    cfg.lpfilter = 'yes';
    cfg.lpfreq = 45;
    cfg.bsfilter = 'yes';
    cfg.bsfreq = [48 52; 96 104];
    cfg.trl = [start stop offset];
    cfg.trl = round(cfg.trl);
    cfg.reref = 'yes';    
    cfg.refchannel = {'M1','M2'};
    preproc = ft_preprocessing(cfg);
    preprocs_ica{NapNumber} = preproc;
end

%APPEND DATA FOR REJECTION AND ICA
for NapNumber=1:numel(preprocs_ica)
cfg = [];
cfg.keepsampleinfo='yes';
preprocs_ica{NapNumber}= ft_appenddata(cfg,preprocs_ica{NapNumber},eog_chans{NapNumber},emg_chans{NapNumber},ecg_chans{NapNumber});
end

%SEGMENT FOR REJECTION AND ICA
for NapNumber = 1:numel(preprocs_ica)
cfg = [];
cfg.length = 5;
preprocs_ica_seg{NapNumber} = ft_redefinetrial(cfg, preprocs_ica{NapNumber});
end

%SEGMENT RAW DATA
for NapNumber = 1:numel(preprocs_raw)
cfg = []
cfg.length = 5
preprocs_raw_seg{NapNumber} = ft_redefinetrial(cfg, preprocs_raw{NapNumber});
end

%SEGMENT NO-FLTR EOG
for NapNumber=1:numel(eog_chans_nofltr)
cfg = [];
cfg.length = 5;
eog_chans_nofltr{NapNumber} = ft_redefinetrial(cfg, eog_chans_nofltr{NapNumber});
end

save ('25LDES_cleaning', 'preprocs_ica','preprocs_ica_seg')
%%
%Reject Artifacts
cfg = [];
cfg.layout = lay;
cfg.method ='trial';
cfg.alim = 100;
data.clean{1,1}=ft_rejectvisual(cfg, preprocs_ica_seg{1,1}); %change numbers for each recording
save('25LDES_noarts','data');

% Remove Artifacts from Raw
cfg = []
preprocs_raw_seg_cl{1,7}=ft_rejectartifact(data_clean7.cfg, preprocs_raw_seg{1,7});


%Remove Artifacts from
%% RUN ICA
for NapNumber=1:numel(data.clean)
cfg =[]
cfg.channel = ({'all', '-V*','-E*' ,'-H*', 'M*'});
data_ica{NapNumber}=ft_selectdata(cfg, data.clean{NapNumber});
end
save ('25LDES_ica','data_ica')

load('P:\3013077.01\Scripts_301307701\Shervin\topography\cuetopo_new.mat')
for NapNumber  = 1:numel(data_ica);
cfg            = [];
cfg.channel = {'all', '-M*'};
cfg.method     = 'runica';
data_ica_comp{NapNumber} = ft_componentanalysis(cfg, data_ica{NapNumber});
end 
save('ICA_res','data_ica_comp');

%% Reject Segments from EOGs
cfg =[]
eog_chans_nofltr{1,7}= ft_rejectartifact(data_clean7.cfg, eog_chans_nofltr{1,7});


%Topoplot ICA
figure;
cfg           = [];
cfg.component = [4 8 15 22 26 28 33 34 35 36 38 45 48 50 51 56 59 65 68 72 73 99 101] ;      % specify the component(s) that should be plotted
cfg.layout    = lay; % specify the layout file that should be used for plotting
cfg.comment   = 'no';
ft_topoplotIC(cfg, data_ica_comp{1,4});

%Databrowse ICA
cfg =[]; 
cfg.viewmode ='component';
cfg.channel = [50 51 56 59 65 68 72 73 99 101]
cfg.zlim = 'maxabs'
cfg.layout = lay;
ft_databrowser(cfg, data_ica_comp{1,4})

cfg= []
cfg.ylim = [-50 50]
cfg.viewmode = 'vertical';
ft_databrowser(cfg, eog_chans_nofltr{1,4});

%% Specify Components to Remove and Remove from Raw Data
rmvComp{1,1} = [4 13 14 18 33 46 51 58 78 92];
rmvComp{1,2} = [19 20 29 34 39 58 68];
rmvComp{1,3} = [4 12 19 24 29 56 74 75 88];
rmvComp{1,4} = [4 8 15 22 56 68 72 101];
rmvComp{1,5} = [2 8 20 31 95 105];
rmvComp{1,6} = [5 7 11 16 30 41 51];
rmvComp{1,7} = [1 6 18 17 42 45 85 105];

for NapNumber=1:numel(preprocs_raw_seg_cl);
cfg =[]
cfg.component = rmvComp{NapNumber}
preprocs_rmvcomps{NapNumber} = ft_rejectcomponent(cfg, data_ica_comp{NapNumber}, preprocs_raw_seg_cl{NapNumber})
end

%% INTERPOLATE BAD CHANNELS

load('P:\3013077.01\Scripts_301307701\Shervin\topography\cue_topo.mat');


for NapNumber=1:numel(preprocs_rmvcomps)
cfg = [];
cfg.method = 'spline';
cfg.trials = 'all';
cfg.badchannel =sub.data.badchannel{NapNumber}';
cfg.layout = lay;
cfg.neighbours = neighbours;
repaired_data{NapNumber} = ft_channelrepair(cfg,preprocs_rmvcomps{NapNumber});
end
save('25LDES_rep', 'repaired_data')
    
%% PREPROCESS THE DATA AND REAPPEND EXT. ELECS.

load('25LDES_rep');
for NapNumber = 1:numel(repaired_data);
    cfg = [];
    cfg.demean = 'yes';
    cfg.detrend = 'yes';
    cfg.reref = 'yes';
    cfg.refchannel    = 'all';
    cfg.refmethod = 'avg';
    cfg.channel = 'EEG';
    cfg.hpfilter = 'yes';
    cfg.lpfilter = 'yes';
    cfg.hpfreq   = 0.3;
    cfg.lpfreq   = 45;
    cfg.bsfilter = 'yes';
    cfg.bsfreq = [48 52; 96 104];
    Data_fltr{NapNumber} = ft_preprocessing(cfg,repaired_data{NapNumber});
end

save('25LDES_data', 'Data_fltr')

%% COMBINE EEG + EXT. ELECS 
%APPEND
for NapNumber=1:numel(Data_fltr)
cfg = [];
ICA_compare{NapNumber}= ft_appenddata(cfg, ICA_compare{NapNumber},eog_chans_nofltr{NapNumber});
end

save ('25LDES_preproc', 'preprocs');

%% REDEFINE TRIALS TO 10s w/ 80% overlap

%% POWER SPECTRAL ANALYSIS

%% STATS

%% FREQ ANALYSIS
cfg.channel = 'all'
cfg.method = 'mtmfft'
cfg.foilim     = [0.3 100];
cfg.tapsmofrq  = 0.003;
preprocs_freq= ft_freqanalysis(cfg, Data_fltr{1,1});
plot(preprocs_freq.freq, preprocs_freq.powspctrm)
xlabel('Frequency (Hz)');
ylabel('absolute power (uV^2)');
legend(preprocs_freq.label);
%% DATA BROWSER
for N=1:numel(Data_fltr_eyes)
cfg = [];
cfg.viewmode = 'vertical';
cfg.ylim = [-12.5 12.5]
cfg.channel = [1:10]
cfg.channelclamped='V*', 'H*';
ft_databrowser(cfg, Data_fltr_eyes{N});
end

cfg = [];
cfg.viewmode = 'vertical';
cfg.ylim = [-25 25]
cfg.channel = [1:10]
cfg.channelclamped = ({'VEOG1', 'HEOG1'});
ft_databrowser(cfg, ICA_compare{1,5});

cfg = [];
cfg.viewmode = 'vertical';
cfg.ylim = [-50 50]
cfg.channel = [1:10]
cfg.channelclamped = ['HEOG1'];
ft_databrowser(cfg, Data_fltr_eyes{1,1});

