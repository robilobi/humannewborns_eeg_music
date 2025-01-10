function  [avgOc1, avgOc2, avgSc1, avgSc2, avgrand, perc_trremoved] = fun_load_epoch_eeg_erp_new(file, allinfo, midi_names,code_corresponding_to_midi, sbj)

load(file);

    load('Layout_Newborn_EEG.mat');
    cfg = [];
    cfg.layout          = lay;
    %cfg.feedback='yes';
    layout              = ft_prepare_layout(cfg);
    cfg                 = [];
    cfg.layout          = layout;
    cfg.method          = 'distance'; % for prepare_neigh
    cfg.neighbourdist   = .25;         % results in avg 5 channels
    cfg.neighbours      = ft_prepare_neighbours(cfg);
    cfg.method          = 'nearest';

    while find(isnan(data_clean.trial{1}(:,:))) >0
        missing = data_clean.label(find(isnan(data_clean.trial{1}(:,1))))';
        cfg.missingchannel  = missing;
        cfg.senstype        = 'eeg';
        cfg.method          = 'average';
        data_clean   = ft_channelrepair(cfg, data_clean);
    end

%% [all_id pitch onsets alliti allioi S0 E0 Sp Ep cond];
%      1     2     3      4       5     6  7  8  9   10
%%%% EPOCHING
baseline_dur = 5; % IN S
new_sampling = data_clean.fsample; % in Hz
n_offsets    = [.100 0.500]; % in s
info = allinfo(:, [1 3 10]); % midi id, onset, condition
[data_out] = Robs_MiniTrialMaker_flagNotes_new(data_clean,midi_names,code_corresponding_to_midi,baseline_dur,new_sampling,n_offsets, info);


idxsubj = size(data_out.trialinfo,2)+1;
data_out.trialinfo(:,idxsubj) = sbj;

%%% Baseline
[ data_out] = Giac_Manual_BaseCorr(data_out, [-50 0]);

%% REMOVE BAD TRIALS
%%% Find bad trials (2SD from mean)
mateeg = nt_trial2mat(data_out.trial);
mateeg(isnan(mateeg)) = 0;
% Reject trials that deviate from the mean by more than twice the average deviation from the mean
good_trials=nt_find_outlier_trials(mateeg, 2); %  nt_find_outlier_trials requires input structure time * channels * trials)
alltrials = 1:size(mateeg, 3);
bad_trials = alltrials(~ismember(alltrials, good_trials));
perc_trremoved = (length(bad_trials)*100)/length(alltrials);

cfg = [];
cfg.trials = good_trials;
data_out = ft_selectdata(cfg, data_out); % amend data to only keep good trials


%%% Select notes high vs low S
cfg =[];
cfg.trials = data_out.trialinfo(:,5) == 1;  % LOW IC
datac1 = ft_selectdata(cfg,data_out);
cfg =[];
cfg.trials = data_out.trialinfo(:,5) == 2;  % HIGH IC
datac2 = ft_selectdata(cfg,data_out);


%%% Divide by real vs shuffled music
cfg =[];
cfg.trials =ismember( datac1.trialinfo(:,2),[1:10]);%   cfg.trials =ismember( datac1.trialinfo(:,1),1:10);
orc1 = ft_selectdata(cfg,datac1);
cfg =[];
cfg.trials =ismember( datac2.trialinfo(:,2),[1:10]);%   cfg.trials =ismember( datac1.trialinfo(:,1),1:10);
orc2 = ft_selectdata(cfg,datac2);
cfg =[];
cfg.trials = ismember(datac1.trialinfo(:,2), [11 15 18 20]);
shc1 = ft_selectdata(cfg,datac1);
cfg =[];
cfg.trials = ismember(datac2.trialinfo(:,2), [11 15 18 20]);
shc2 = ft_selectdata(cfg,datac2);

%%% Create surrogate ERP (random trigger)
baseline_dur = 5; % IN S
new_sampling = data_clean.fsample; % in Hz
n_offsets    = [.100 .500]; % in s
[data_out] = Robs_MiniTrialMaker_random_new(data_clean,midi_names,code_corresponding_to_midi,baseline_dur,new_sampling,n_offsets, info);
cfg =[];
cfg.trials = randperm( length(data_out.trial), 3000);
datarandom = ft_selectdata(cfg,data_out);
datarandom = Giac_Manual_BaseCorr( datarandom, [-50 0]);


%%%% Average over trials 
cfg = [];
avgOc1 = ft_timelockanalysis(cfg, orc1);
avgOc2 = ft_timelockanalysis(cfg, orc2);
avgSc1 = ft_timelockanalysis(cfg, shc1);
avgSc2 = ft_timelockanalysis(cfg, shc2);
avgrand = ft_timelockanalysis(cfg, datarandom);


end

