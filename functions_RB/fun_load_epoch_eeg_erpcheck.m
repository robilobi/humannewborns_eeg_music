function  [avgo, avgs] = fun_load_epoch_eeg_erpcheck(file, midi_names,code_corresponding_to_midi, midi_path, sbj)

load(file);

% allinfo [all_id allpitches allonsets alliti allioi S Sp So Cond ];
% allinfo    1       2         3         4     5     6  7  8   9

%%%% EPOCHING
baseline_dur = 5; % IN S
new_sampling = data_clean.fsample; % in Hz
n_offsets    = [.500 1.500]; % in s
[data_out] = RFT_MiniTrialMaker(data_clean,midi_names,code_corresponding_to_midi,midi_path,baseline_dur,new_sampling,n_offsets);

clear data_int
%     [ nan_trials ] = Giac_findNanTrials( data_out, 'OnlyTrial' ); % session 12 contains some NaN
%     [data_out] = Giac_removeTrials(data_out,nan_trials,'reject');
idxsubj = size(data_out.trialinfo,2)+1;
data_out.trialinfo(:,idxsubj) = sbj;
% data_out.trialinfo = melNumber, melID, Session, Note Number, Surprise Condition, IC, ITI, IOI, subject


%%% Baseline
[ data_out] = Giac_Manual_BaseCorr( data_out, [-100 0]);

%%% Divide by real vs shuffled music
cfg =[];
cfg.trials =ismember( data_out.trialinfo(:,2),[1:10]);%   cfg.trials =ismember( datac1.trialinfo(:,1),1:10);
or = ft_selectdata(cfg,data_out);

cfg =[];
cfg.trials = ismember(data_out.trialinfo(:,2), [11 15 18 20]);
sh = ft_selectdata(cfg,data_out);


cfg = [];
avgo = ft_timelockanalysis(cfg, or);
avgs = ft_timelockanalysis(cfg, sh);


end

