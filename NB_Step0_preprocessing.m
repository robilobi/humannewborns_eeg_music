%% Step 1 Bach Newborns: Load EEG and preprocess
% Steps are Bandpass filtering, Bad channel rejection, ASR,
% Automatic ICA, Bad channel rejection part 2 and interpolation, 
% Rereferencing 
%%%%%%%%%%%%%% Roberta Bianco Nov 2023

%% initialise fieldtrip
% clc
% clear all
% close all

U = 'C:\Users\robianco\OneDrive - Fondazione Istituto Italiano Tecnologia';
addpath([U '\MATools\fieldtrip-20220321\']); ft_defaults;
addpath([U '\MATools\eeglab2022.1\' ]);
eeglab
addpath('.\functions_RB')
addpath('.\templates')
ft_defaults;

dataINFolder = ([U '\BACHNB\DATA\']);
dataOUTfolder  =  ([U '\BACHNB\RESULTS_NB\']); mkdir(dataOUTfolder);
wav_path = ([U '\BACHNB\STIMULI\test_wav_final_beep\']);
eegsubjects = dir([dataINFolder,'s*']);

%% get end of trials from wav files
stim_names = {'audio01','audio02','audio03','audio04','audio05','audio06','audio07','audio08',...
    'audio09','audio10', 'shf01','shf05','shf08','shf10'};
durations = RB_getsonglength(wav_path, stim_names); % in seconds

%% start loop
subjidx = [1:25 27:41 43:50 55:64]; %those not included have no events FINAL N = 58
nameoutfile = '1-30Hz';

for k = subjidx

    % load names of all files
    cd([dataINFolder eegsubjects(k).name])
    FileNames      = {dir(['*' '.vhdr']).name} ;

    % get baby ID
    baby_name      = eegsubjects(k).name;

    % load auditory conditions
    triggers = {'S  1','S  2','S  3','S  4','S  5','S  6','S  7','S  8','S  9', 'S 10','S 11','S 15','S 18','S 20'};
    trg = [1:10 11 15 18 20]; % map triggers to integers
    data_all  = RB_loadConditionData(FileNames{1},0, durations, triggers, trg);

    % bandpass filter data
    data_filt      = RB_filterEEG(data_all,'bp',[1 30]); %BETA (13-30Hz), ALPHA (8-12 Hz)

    % remove trials that were not recorded
    [ nan_trials ] = Giac_findNanTrials( data_filt, 'OnlyTrial' );
    [data_rmtr] = Giac_removeTrials(data_filt,nan_trials,'reject');

    %RFT_EEG_Visualise( data_rmtr, 'all', [-30 30], 'Layout_Newborn_EEG.mat')

    %% Exclude bad channels and trials
    % check for flatlines and noisy channels according to SD, mean and peak to peak
    steps_to_do    = {'flat','corr-noise'};
    criteria       = {[5],[0.1,20]};
    flat_channels  = RFT_detect_FlatCorrNoise_electrodes(data_rmtr,'Layout_Newborn_eeglab.mat',steps_to_do,criteria);
    data_flat      = Giac_removeChannels(data_rmtr, flat_channels.flat );
    data_flat      = Giac_removeChannels(data_flat, flat_channels.corrNoise );

    % get noisy bad channels according to SD,mean and peak to peak info
    bad_channels   = Giac_EEG_CatchNoisyElectrodes(data_flat, 'all', 2.75, 'recursive' );
    data_badchan   = Giac_removeChannels(data_flat, bad_channels );

    %% reref
    cfg                      = [];
    cfg.channel              = {'EEG'};
    cfg.reref                = 'yes';
    cfg.refchannel           = {'F9', 'P9', 'Iz', 'P10', 'F10'}; 
    data_reref               = ft_preprocessing(cfg,data_badchan);

    %% Artefact correction
     data_asr       = RFT_clean_asr_combined_trials_ftstruct(data_reref,'Layout_Newborn_eeglab.mat',  5,[],[],[]);

    %% ICA
    [data_ICA,rejected_comps] = RFT_IClabel(data_asr,'Layout_Newborn_eeglab.mat',30,[0 0;0 0; 0.5 1; 0 0; 0 0; 0 0; 0 0]);
    close all

    %% DOWNSAMPLE AT 100 HZ (ENOUGH FOR ERP AND TRF ANALYSIS)
    cfg                      = [];
    cfg.resamplefs           = 100;
    cfg.detrend              = 'no';
    cfg.demean               = 'no';
    [data_dw]             = ft_resampledata(cfg, data_ICA);

    %% REPAIR CHANNELS
    bad_channels2     = Giac_EEG_CatchNoisyElectrodes(data_dw, 'all', 2.75, 'recursive' );
    bad_channels      = [bad_channels bad_channels2];
    data_prepint      = Giac_removeChannels( data_dw, bad_channels );

    cfg                 = [];
    cfg.layout          = layout;
    cfg.method          = 'distance'; % for prepare_neigh
    cfg.neighbourdist   = 0.18;         % results in min 4 channels, mean 8.6
    cfg.neighbours      = ft_prepare_neighbours(cfg);     
    cfg.method          = 'nearest';
    cfg.senstype        = 'eeg';
    %cfg.elec            = ft_read_sens('standard_1020.elc'); % 
    cfg.badchannel      = bad_channels;
    data_clean          = ft_channelrepair(cfg, data_prepint);

%% if there are still missing channels after interpolation redo it
    cfg                 = [];
    cfg.layout          = layout;
    cfg.method          = 'distance'; % for prepare_neigh
    cfg.neighbourdist   = .18;         % results in min 4 neighbours,  mean 8.6 
    cfg.neighbours      = ft_prepare_neighbours(cfg);

    while find(isnan(data_clean.trial{1}(:,:))) >0
        missing = data_clean.label(find(isnan(data_clean.trial{1}(:,1))))';
        cfg.missingchannel  = missing;
        cfg.senstype        = 'eeg';
        cfg.method          = 'spline';
        data_clean   = ft_channelrepair(cfg, data_clean);
    end

    %%% STORE removed stuff
    Bad_channels{k} =  [bad_channels(:)'];
    Bad_IcaCmp{k}   = rejected_comps;

    %% save data
    save(strcat(dataOUTfolder,baby_name,'_PROC_', nameoutfile),'data_clean');
end
save(strcat(dataOUTfolder,'allbabies','_META'),'Bad_channels', 'Bad_IcaCmp');


