function [EEG] = convertft2eeglab(ft_struct,eeglab_template_file)
load(ft_struct);
ft_struct = data_clean;

all_trials = cat(2,ft_struct.trial{:});
eeglab_template = load(eeglab_template_file);
EEG.data = all_trials;

%%% update other fields of interest %%%
[EEG.nbchan , EEG.pnts] = size(all_trials);  % number of channels , number of time points
EEG.trials = size(ft_struct.trial,2);                                        % because 1 big trial (of concatenated data)
EEG.srate = ft_struct.fsample;                         % important if you have a sampling rate ≠ 2048Hz
EEG.chanlocs = eeglab_template.eeglab_template.chanlocs;
EEG.ref = 'Cz';

%%% update chanlocs if you have removed channels %%%
chans_removed_in_ft=[];
for ch=1:length(EEG.chanlocs)
    if ~ismember(EEG.chanlocs(ch).labels,ft_struct.label)   
        chans_removed_in_ft(end+1)=ch;        end
end
EEG.chanlocs(chans_removed_in_ft) = [];

%%% reorder chanlocs if the channels are ordered differently in ft struct
%%% (can happen for ex if you interpolated channels)
assert(length(ft_struct.label)==length(EEG.chanlocs),'eeglab chanlocs should have same number of channels that ft struct');
for ch=1:length(ft_struct.label)
    new_idx(ch) = find(strcmp(ft_struct.label(ch) , {EEG.chanlocs.labels}));
end
EEG.chanlocs = EEG.chanlocs(new_idx);

end