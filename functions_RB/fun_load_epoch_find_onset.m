function  [avg, avgrand] = fun_load_epoch_find_onset(file, midi_names,code_corresponding_to_midi, midi_path, sbj)

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


%%%% EPOCHING
baseline_dur = 5; % IN S
new_sampling = data_clean.fsample; % in Hz
n_offsets    = [.050 0.500]; % in s
data = RFT_MiniTrialMaker(data_clean,midi_names,code_corresponding_to_midi,midi_path,baseline_dur,new_sampling,n_offsets);
idxsubj = size(data.trialinfo,2)+1;
data.trialinfo(:,idxsubj) = sbj;
data = Giac_Manual_BaseCorr(data, [-50 0]);

%%% Create surrogate ERP (random trigger)
baseline_dur = 5; % IN S
new_sampling = data_clean.fsample; % in Hz
n_offsets    = [.050 0.500]; % in s
datarandom = RFT_MiniTrialMaker_rnd(data_clean,midi_names,code_corresponding_to_midi,midi_path,baseline_dur,new_sampling,n_offsets);
datarandom.trialinfo(:,idxsubj) = sbj;
datarandom = Giac_Manual_BaseCorr( datarandom, [-50 0]);

%%%% Average over trials
cfg = [];
avg = ft_timelockanalysis(cfg, data);
avgrand = ft_timelockanalysis(cfg, datarandom);

end

