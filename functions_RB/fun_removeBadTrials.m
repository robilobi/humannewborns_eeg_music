function [data, good_trials] = fun_removeBadTrials(data, threshold)
%%
% Remove trials exceeding 2SD from mean
%% Roberta Bianco
time = data.time{1};
x = cat(3,data.trial{:});
x(isnan(x)) = 0;
x = permute(x,[2 1 3]); % nt_find_outlier_trials requires input structure time * channels * trials)
xb = nt_demean(x, find(time<0)); %baseline
 % Reject trials that deviate from the mean by more than twice the average deviation from the mean
good_trials=nt_find_outlier_trials(xb, threshold); % creates array with good trials based on threshold
cfg = [];
cfg.trials = good_trials;
data = ft_selectdata(cfg, data); % amend data to only keep good trials

end

