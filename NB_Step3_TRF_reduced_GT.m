% Run mTRF analysis on the EEG data from one subject
% MODELS to run: Full, Full with random music predictors, So, Sp, ioi, iti
%%%%%%%%%%%%%%%%% Roberta Bianco, Nov 2024, Rome %%%%%%%%%%%%%%%
clear all
close all
%%% PATHS
U = 'C:\Users\robianco\OneDrive - Fondazione Istituto Italiano Tecnologia';
addpath([U '\MATools\fieldtrip-20220321\']); ft_defaults;
addpath('.\functions_RB');
addpath('.\templates');
ft_defaults;
addpath([U '\MATools\mTRF-Toolbox-master\mtrf']); % add the mTRF toolbox functions

stim_path = ([U '\BACHNB\STIMULI\']);
dataINFolder = ([U '\BACHNB\RESULTS_NB\']);
dataOUTfolder = [U '\BACHNB\RESULTS_NB\TRF_GT_stm_rnd_sosp\']; mkdir(dataOUTfolder);
matrixtoload ='FINAL_StimMatrix_both+';


stim_delay = 5; % delay of the actual stimulus start relative to the start of the EEG trial
skip_time = 5.2; % amount of time to exclude from the beginning and ending of the stimulus and EEG before modeling

% Model parameters
Fs = 100;
tmin = -50;
tmax = 400;
lambdas = [0 10.^(-4:8)];

%%% SUBJECT POOL
eegsubjects = dir([dataINFolder,'*_PROC_1-30Hz*']);
subjidx = 1:length(eegsubjects);
subjidx = [1:10 12:33 35:44 46:47 54:58]; % those excluded don't have all melodies


stimnames = {'audio01','audio02','audio03','audio04','audio05', ...
    'audio06','audio07','audio08', 'audio09','audio10',...
    'shf01','shf05','shf08','shf10'};
stimindx = 1:14;

cfg = [];
cfg.layout = 'Layout_Newborn_EEG.mat';
cfg.baselinetype = 'absolute';
cfg.skipscale   = 'yes'; cfg.skipcomnt   = 'yes' ;
layout = ft_prepare_layout(cfg);
norm_chanlbls =layout.label([2:26 28:65]);


%% start
for sub = subjidx

    %%% Load the EEG data
    file = [dataINFolder eegsubjects(sub).name];
    disp('Loading EEG data...');
    load(file);

    %%% 
    alleeg = {}; stim_vals = [];
    for t = 1:length(data_clean.trial)
        alleeg = [alleeg; data_clean.trial{t}'];
    end
    stim_vals = data_clean.trialinfo(:,1);
    Fs = data_clean.fsample;
    chan_lbls = data_clean.label;

    % Get the stimulus names based on the trigger values
    allstim_names = get_stims_from_triggers(stim_vals);
    stim_played = unique(allstim_names);


    % Average across repetitions of melodies
    [eegavg]= fun_avgtrialsbyname(alleeg,stim_played,allstim_names);

    % Reorder channels according to template channels
    [a,ch] = ismember(norm_chanlbls, chan_lbls);
    for m = 1:length(stim_played) % sub 11 lacks song 9
        for c = 1:length(ch)
            neweegavg{m}(:,c) = eegavg{m}(:, ch(c));
        end
    end
    eegavg = neweegavg;

    % Normalise eeg by the global variance over time (concat trials)
    EEG_allTrials_concat = vertcat(eegavg{:}); % concat trials x chan
    EEG_global_std = sqrt(mean(var(EEG_allTrials_concat,[],1)));

    % Load the stim matrix
    load([stim_path matrixtoload]);

    % Some songs are missing ion some babies (e.g., sub 11, 34). Remove
    % songs that were not played, and are not in the eeg from the stim matrix
    if length(stim_played) < length(stimnames)
        mismatches = find(ismember(stimnames, stim_played)==0);
        for s = mismatches; stim_all(s) = [];end
    end

    % load groundtruth data (GAVG data)
    load([ 'GroundTruthEEG_NB.mat']);
    % Remove amount of time from the start and end of the stimulus and EEG
    if skip_time > 0
        fprintf('Removing %.1f s from start and end of the stimulus and EEG...\n',skip_time);
        for s =1:length(stim_played)
            nidx = size(stim_all{s},1); % number of time indexes in the stimulus
            use_idx = (skip_time*Fs+1):(nidx - 4*Fs);
            stim_all{s} = stim_all{s}(use_idx,:);
            d = detrend(eegavg{s}(use_idx,:),0); % also shift the mean to 0
            eegavg{s} = (d ./ EEG_global_std);
            dd = meanmegasubjeeg{s}(use_idx,:); % also shift the mean to 0
            meanmegasubjeegcut{s} = dd;
            % normalize the stimulus inputs so their rms = 1
            stim_all{s} = stim_all{s}./(ones(length(use_idx),1)*rms(stim_all{s}));
        end
    end


    %% Quantify prediction accuracy 
    %%% (20-4-2022) Currently this is training on all stimuli presented
    pred_names = {'Spflux', 'onset', 'ITI', 'IOI', ...
    'Soioiratio', 'Eoioiratio', ...
    'Sppitch','Eppitch' };

    %% TRF - Compute model weights
    nreg = length(pred_names);
    allregr = 1:nreg;
    k=1;
    % full model
    keepregr = [1 2 3 4 5 6 7 8]; %so sp
%     keepregr = [1 2 3 4 11 12 9 10]; %soc spc
    regressors = stim_all;for j = 1:length(stim_all); regressors{j}(:,setdiff(allregr, keepregr)) = [];end
    [stats,m] = itercvtest_GT(regressors,eegavg,Fs,1,tmin,tmax,lambdas, meanmegasubjeegcut);
    stats.stim_names = stim_played;
    model(k).statsall{sub} =stats;
    model(k).mdl{sub} = m;
    model(k).predname = pred_names(keepregr);
    k = k+1;
    
    % random MUSIC / acoustics
    rndregr  = [0 0 0 0 1 1 1 1]; % regressors to randomise
    regressors = stim_all;for j = 1:length(stim_all); regressors{j}(:,setdiff(allregr, keepregr)) = [];end
    for h = 1:length(regressors) % for each melody
        variables = regressors{h};
        [var_sh] = fun_TRFmake_idyom_random(variables, find(rndregr==1));
        regressors{h} =var_sh;
    end
    [stats,m] = itercvtest_GT(regressors,eegavg,Fs,1,tmin,tmax,lambdas, meanmegasubjeegcut);
    stats.stim_names = stim_played;
    model(k).statsall{sub} =stats;
    model(k).mdl{sub} = m;
    model(k).predname = pred_names(keepregr);
    k = k+1;

    % random So
    rndregr  = [0 0 0 0 1 1 0 0]; % regressors to randomise
    regressors = stim_all;for j = 1:length(stim_all); regressors{j}(:,setdiff(allregr, keepregr)) = [];end
    for h = 1:length(regressors) % for each melody
        variables = regressors{h};
        [var_sh] = fun_TRFmake_idyom_random(variables, find(rndregr==1));
        regressors{h} =var_sh;
    end
    [stats,m] = itercvtest_GT(regressors,eegavg,Fs,1,tmin,tmax,lambdas, meanmegasubjeegcut);
    stats.stim_names = stim_played;
    model(k).statsall{sub} =stats;
    model(k).mdl{sub} = m;
    model(k).predname = pred_names(keepregr);
    k = k+1;

    % random Sp
    rndregr  = [0 0 0 0 0 0 1 1]; % regressors to randomise
    regressors = stim_all;for j = 1:length(stim_all); regressors{j}(:,setdiff(allregr, keepregr)) = [];end
    for h = 1:length(regressors) % for each melody
        variables = regressors{h};
        [var_sh] = fun_TRFmake_idyom_random(variables, find(rndregr==1));
        regressors{h} =var_sh;
    end
    [stats,m] = itercvtest_GT(regressors,eegavg,Fs,1,tmin,tmax,lambdas, meanmegasubjeegcut);
    stats.stim_names = stim_played;
    model(k).statsall{sub} =stats;
    model(k).mdl{sub} = m;
    model(k).predname = pred_names(keepregr);
    k = k+1;

    % random IOI
    rndregr  = [0 0 0 1 0 0 0 0]; % regressors to randomise
    regressors = stim_all;for j = 1:length(stim_all); regressors{j}(:,setdiff(allregr, keepregr)) = [];end
    for h = 1:length(regressors) % for each melody
        variables = regressors{h};
        [var_sh] = fun_TRFmake_idyom_random(variables, find(rndregr==1));
        regressors{h} =var_sh;
    end
    [stats,m] = itercvtest_GT(regressors,eegavg,Fs,1,tmin,tmax,lambdas, meanmegasubjeegcut);
    stats.stim_names = stim_played;
    model(k).statsall{sub} =stats;
    model(k).mdl{sub} = m;
    model(k).predname = pred_names(keepregr);
    k = k+1;

    % random ITI
    rndregr  = [0 0 1 0 0 0  0 0]; % regressors to randomise
    regressors = stim_all;for j = 1:length(stim_all); regressors{j}(:,setdiff(allregr, keepregr)) = [];end
    for h = 1:length(regressors) % for each melody
        variables = regressors{h};
        [var_sh] = fun_TRFmake_idyom_random(variables, find(rndregr==1));
        regressors{h} =var_sh;
    end
    [stats,m] = itercvtest_GT(regressors,eegavg,Fs,1,tmin,tmax,lambdas, meanmegasubjeegcut);
    stats.stim_names = stim_played;
    model(k).statsall{sub} =stats;
    model(k).mdl{sub} = m;
    model(k).predname = pred_names(keepregr);
    clear eegavg neweegavg

end

fn= 'TRFout_newborns_AM.mat';
save([dataOUTfolder fn], 'model', 'stim_all',...
    'stimnames', 'norm_chanlbls', '-mat');

