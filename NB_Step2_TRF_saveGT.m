% Create a ground truth eeg dataset from the average of all participats. It
% will be used in the TRF analysis
%%%%%%%%%%%%%%% Roberta Bianco Nov 2024, Rome %%%%%%%%%%%%%%%%%%%%%

clear all
close all
%%% PATHS
U = 'C:\Users\robianco\OneDrive - Fondazione Istituto Italiano Tecnologia';
addpath([U '\MATools\fieldtrip-20220321\']); ft_defaults;
addpath('.\functions_RB')
addpath('.\templates')

dataINFolder = ([U '\BACHNB\RESULTS\']);
dataOUTfolder = [U '\BACHNB\RESULTS\TRF_bysubj_1-30Hz_reduced2_GT\']; mkdir(dataOUTfolder);

%%% SUBJECT POOL
eegsubjects = dir([dataINFolder,'*_PROC_1-30Hz*']);
subjidx = [1:10 12:33 35:44 46:47 54:58];

stimnames = {'audio01','audio02','audio03','audio04','audio05', ...
    'audio06','audio07','audio08', 'audio09','audio10',...
    'shf01','shf05','shf08','shf10'};

cfg = [];
cfg.layout = 'Layout_Newborn_EEG.mat';
cfg.baselinetype = 'absolute';
cfg.skipscale   = 'yes'; cfg.skipcomnt   = 'yes' ;
layout = ft_prepare_layout(cfg);
norm_chanlbls =layout.label([2:26 28:65]);


%%% START
megasubjeeg = cell(length(stimnames), 1);
countsub = 1;
for sub = subjidx

    % Load the EEG data
    file = [dataINFolder eegsubjects(sub).name];
    disp('Loading EEG data...');
    load(file);

    % Organise in TRF format
    alleeg = {}; stim_vals = [];
    for t = 1:length(data_clean.trial)
        alleeg = [alleeg; data_clean.trial{t}'];
    end
    stim_vals = data_clean.trialinfo(:,1);
    chan_lbls = data_clean.label;

    % Get the stimulus names based on the trigger values
    allstim_names = get_stims_from_triggers(stim_vals);
    stim_played = unique(allstim_names);

    % Average across repetitions (if any) of melodies
    [eegavg]= fun_avgtrialsbyname(alleeg,stim_played,allstim_names);

    % Normalise eeg by the global variance over time (concat trials)
    EEG_allTrials_concat = vertcat(eegavg{:}); % concat trials x chan
    EEG_global_std = sqrt(mean(var(EEG_allTrials_concat,[],1)));

    % Detrend and normalise
    detrendFunc = @(x) detrend(x, 0);% 
    detrendedeeg = cellfun(detrendFunc, eegavg, 'UniformOutput', false);
    normFunc = @(x) (x ./ EEG_global_std);% 
    normeeg = cellfun(normFunc, detrendedeeg, 'UniformOutput', false);

    % Reorder channels according to template channels
    [a,ch] = ismember(norm_chanlbls, chan_lbls);
    [a, sidx] = ismember(stim_played, stimnames);
    count = 1;
    for m = sidx' % sub 11 lacks song 9
        for c = 1:length(ch)
            megasubjeeg{m}(:,c, countsub) = normeeg{count}(:, ch(c));
        end
        count = count +1;
    end
    countsub = countsub +1;
end

% Average across subjects by melody
meanFunc = @(x) nanmean(x, 3);% Function to calculate the mean along the third dimension for each cell
meanmegasubjeeg = cellfun(meanFunc, megasubjeeg, 'UniformOutput', false);

save([dataOUTfolder 'GroundTruthEEG_NB.mat'], "meanmegasubjeeg", '-mat');
